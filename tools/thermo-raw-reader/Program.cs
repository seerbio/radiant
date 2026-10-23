using System.Collections.Concurrent;
using System.Globalization;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoRawFileParser.Writer;

// This process reads RAW directly. Standard output carries only the binary scan
// protocol; all diagnostics go to standard error. It never writes mzML or XML.
internal static class Program
{
    private const uint ProtocolVersion = 4;
    private const uint UnknownSpectrumCount = uint.MaxValue;
    private const int MaxPointsPerScan = 10_000_000;
    private const int ReaderWorkerLimit = 4;
    private const int MaxScansPerChunk = 64;
    private const int MaxBytesPerChunk = 4 * 1024 * 1024;
    private static readonly byte[] Magic = Encoding.ASCII.GetBytes("RDRRAW01");
    private static readonly byte[] Footer = Encoding.ASCII.GetBytes("RAWDONE1");

    private sealed record ScanRange(long First, long EndExclusive);
    private sealed record SerializedChunk(byte[] Buffer, int Length, int SpectrumCount);

    private sealed class ChunkBuffer : IDisposable
    {
        private MemoryStream _buffer;
        private BinaryWriter _writer;

        public ChunkBuffer()
        {
            _buffer = new MemoryStream();
            _writer = new BinaryWriter(_buffer, Encoding.UTF8, leaveOpen: true);
        }

        public BinaryWriter Writer => _writer;
        public int SpectrumCount { get; private set; }

        public void AddSpectrum() => ++SpectrumCount;

        public bool IsFull => SpectrumCount >= MaxScansPerChunk || _buffer.Length >= MaxBytesPerChunk;

        public void Publish(BlockingCollection<SerializedChunk> chunks)
        {
            _writer.Flush();
            chunks.Add(new SerializedChunk(
                _buffer.GetBuffer(),
                checked((int)_buffer.Length),
                SpectrumCount));
            _writer.Dispose();
            _buffer.Dispose();
            _buffer = new MemoryStream();
            _writer = new BinaryWriter(_buffer, Encoding.UTF8, leaveOpen: true);
            SpectrumCount = 0;
        }

        public void Dispose()
        {
            _writer.Dispose();
            _buffer.Dispose();
        }
    }

    private static int Main(string[] args)
    {
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
        CultureInfo.CurrentUICulture = CultureInfo.InvariantCulture;
        if (args.Length == 1 && args[0] == "--version")
        {
            Console.WriteLine("RadiantThermoRawReader 1.0.0; binary protocol 4");
            Console.WriteLine("ThermoRawFileParser scan utilities and ThermoFisher CommonCore RawFileReader");
            Console.WriteLine("RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved.");
            return 0;
        }
        if (args.Length != 1 || args[0] is "--help" or "-h")
        {
            Console.Error.WriteLine("Usage: RadiantThermoRawReader INPUT.raw");
            Console.Error.WriteLine("Streams centroided MS1/MS2 scans to stdout using Radiant's binary protocol.");
            return args.Length == 1 ? 0 : 2;
        }
        try
        {
            ReadRaw(Path.GetFullPath(args[0]));
            return 0;
        }
        catch (Exception error)
        {
            Console.Error.WriteLine($"Thermo RAW read failed: {error.Message}");
            return 1;
        }
    }

    private static void ReadRaw(string path)
    {
        if (!File.Exists(path))
            throw new FileNotFoundException("RAW file does not exist", path);

        using var primaryRaw = OpenRaw(path);
        long firstSpectrum = primaryRaw.RunHeaderEx.FirstSpectrum;
        long lastSpectrum = primaryRaw.RunHeaderEx.LastSpectrum;
        if (firstSpectrum <= 0 || lastSpectrum < firstSpectrum)
            throw new IOException("The RAW file has an invalid native scan range.");

        long nativeSpectrumCount = lastSpectrum - firstSpectrum + 1;
        int workerCount = (int)Math.Min(ReaderWorkerLimit, nativeSpectrumCount);
        var ranges = new ScanRange[workerCount];
        for (int worker = 0; worker < workerCount; ++worker)
        {
            long rangeStart = firstSpectrum + nativeSpectrumCount * worker / workerCount;
            long rangeEnd = firstSpectrum + nativeSpectrumCount * (worker + 1) / workerCount;
            ranges[worker] = new ScanRange(rangeStart, rangeEnd);
        }

        using var chunks = new BlockingCollection<SerializedChunk>(workerCount * 2);
        using var output = new BufferedStream(Console.OpenStandardOutput(), 1024 * 1024);
        using var writer = new BinaryWriter(output, Encoding.UTF8, leaveOpen: true);
        writer.Write(Magic);
        writer.Write(ProtocolVersion);
        writer.Write(UnknownSpectrumCount);

        var workers = new Task[workerCount];
        for (int worker = 0; worker < workerCount; ++worker)
        {
            ScanRange range = ranges[worker];
            IRawDataPlus? existingReader = worker == 0 ? primaryRaw : null;
            workers[worker] = Task.Factory.StartNew(
                () => ReadRange(path, range, existingReader, chunks),
                CancellationToken.None,
                TaskCreationOptions.LongRunning,
                TaskScheduler.Default);
        }

        Task completion = Task.Run(async () =>
        {
            try
            {
                await Task.WhenAll(workers).ConfigureAwait(false);
            }
            finally
            {
                chunks.CompleteAdding();
            }
        });

        uint spectrumCount = 0;
        foreach (SerializedChunk chunk in chunks.GetConsumingEnumerable())
        {
            writer.Write(checked((uint)chunk.Length));
            writer.Write(checked((uint)chunk.SpectrumCount));
            writer.Write(chunk.Buffer, 0, chunk.Length);
            spectrumCount = checked(spectrumCount + (uint)chunk.SpectrumCount);
        }
        completion.GetAwaiter().GetResult();

        if (spectrumCount == 0)
            throw new IOException("The RAW file contains no MS1 or MS2 spectra.");

        writer.Write(Footer);
        writer.Write(spectrumCount);
        writer.Flush();
    }

    private static IRawDataPlus OpenRaw(string path)
    {
        IRawDataPlus raw = RawFileReaderFactory.ReadFile(path);
        try
        {
            if (!raw.IsOpen)
                throw new IOException("The Thermo library could not open the RAW file.");
            if (raw.IsError)
                throw new IOException($"The RAW file has a Thermo reader error: {raw.FileError}.");
            if (raw.InAcquisition)
                throw new IOException("The RAW file is still being acquired.");
            raw.SelectInstrument(Device.MS, 1);
            // ThermoRawFileParser includes reference/exception centroid data by
            // default. RawFileReader's default differs and can omit valid peaks.
            raw.IncludeReferenceAndExceptionData = true;
            return raw;
        }
        catch
        {
            raw.Dispose();
            throw;
        }
    }

    private static void ReadRange(
        string path,
        ScanRange range,
        IRawDataPlus? existingReader,
        BlockingCollection<SerializedChunk> chunks)
    {
        IRawDataPlus raw = existingReader ?? OpenRaw(path);
        bool ownsReader = existingReader is null;
        try
        {
            using var chunk = new ChunkBuffer();
            for (long currentNumber = range.First; currentNumber < range.EndExclusive; ++currentNumber)
            {
                int nativeNumber = checked((int)currentNumber);
                var scanEvent = raw.GetScanEventForScanNumber(nativeNumber);
                int level = (int)scanEvent.MSOrder;
                if (level is not (1 or 2))
                    continue;

                var (masses, intensities) = ReadCentroids(raw, scanEvent, nativeNumber);
                if (masses.Length != intensities.Length || masses.Length > MaxPointsPerScan)
                    throw new IOException($"Invalid centroid array sizes at scan {nativeNumber}.");

                float precursor = -1, lower = -1, upper = -1, collision = -1;
                if (level == 2)
                {
                    // These are the parser's public scan utilities. Isolation
                    // offsets follow MzMlSpectrumWriter.ConstructPrecursorList.
                    IReaction? reaction = SpectrumWriter.GetReaction(scanEvent, nativeNumber);
                    if (reaction == null)
                        throw new IOException($"MS2 scan {nativeNumber} has no activation reaction.");
                    var trailer = new ScanTrailer(raw.GetTrailerExtraInformation(nativeNumber));
                    double? width = reaction.IsolationWidth >= 0
                        ? reaction.IsolationWidth
                        : trailer.AsDouble("MS2 Isolation Width:");
                    precursor = (float)reaction.PrecursorMass;
                    if (width.HasValue && width.Value >= 0)
                    {
                        double offset = width.Value / 2 + reaction.IsolationWidthOffset;
                        lower = (float)(width.Value - offset);
                        upper = (float)offset;
                    }
                    if (reaction.CollisionEnergyValid)
                        collision = (float)reaction.CollisionEnergy;
                }

                BinaryWriter frameWriter = chunk.Writer;
                frameWriter.Write(checked(36u + (uint)masses.Length * 8u));
                // Protocol 4 keys each frame by its native scan number. Chunks
                // from independent readers can therefore be serialized as
                // soon as they are ready; C++ restores sequential scan numbers.
                frameWriter.Write(nativeNumber);
                frameWriter.Write(nativeNumber);
                frameWriter.Write(level);
                // The mapped mzML reader first parses minutes to float and then
                // multiplies by 60. Preserve that order of operations in C++.
                frameWriter.Write((float)raw.RetentionTimeFromScanNumber(nativeNumber));
                frameWriter.Write(collision);
                frameWriter.Write(precursor);
                frameWriter.Write(lower);
                frameWriter.Write(upper);
                frameWriter.Write(checked((uint)masses.Length));
                for (int point = 0; point < masses.Length; ++point)
                {
                    frameWriter.Write((float)masses[point]);
                    frameWriter.Write((float)intensities[point]);
                }
                chunk.AddSpectrum();

                if (chunk.IsFull)
                {
                    chunk.Publish(chunks);
                }
            }

            if (chunk.SpectrumCount > 0)
                chunk.Publish(chunks);
        }
        finally
        {
            if (ownsReader)
                raw.Dispose();
        }
    }

    private static (double[] Masses, double[] Intensities) ReadCentroids(
        IRawDataPlus raw, IScanEvent scanEvent, int scanNumber)
    {
        var scanStatistics = raw.GetScanStatsForScanNumber(scanNumber);
        if (scanStatistics.IsCentroidScan
            && scanStatistics.SpectrumPacketType == SpectrumPacketType.FtCentroid)
        {
            // Read native FT centroid data directly instead of constructing a
            // Scan object for every spectrum. This is the RawFileReader fast
            // path and returns the same centroid stream selected by Scan.FromFile.
            var centroidStream = raw.GetCentroidStream(scanNumber, true);
            return (centroidStream.Masses, centroidStream.Intensities);
        }

        // Preserve the parser's scan-object path for other packet types,
        // including segmented data and profile scans that need centroiding.
        var scan = Scan.FromFile(raw, scanNumber);
        double[] masses;
        double[] intensities;
        if (scan.HasCentroidStream)
        {
            masses = scan.CentroidScan.Masses;
            intensities = scan.CentroidScan.Intensities;
        }
        else if (scan.SegmentedScan.PositionCount > 0)
        {
            var segmented = scanEvent.ScanData == ScanDataType.Profile
                ? Scan.ToCentroid(scan).SegmentedScan
                : scan.SegmentedScan;
            masses = segmented.Positions;
            intensities = segmented.Intensities;
        }
        else
        {
            return (Array.Empty<double>(), Array.Empty<double>());
        }
        return (masses, intensities);
    }
}
