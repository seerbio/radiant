using System.Globalization;
using System.Text;
using ThermoFisher.CommonCore.Data;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoRawFileParser.Writer;

// This process reads RAW directly. Standard output carries only the binary scan
// protocol; all diagnostics go to standard error. It never writes mzML or XML.
internal static class Program
{
    private const uint ProtocolVersion = 3;
    private const uint UnknownSpectrumCount = uint.MaxValue;
    private const int MaxPointsPerScan = 10_000_000;
    private static readonly byte[] Magic = Encoding.ASCII.GetBytes("RDRRAW01");
    private static readonly byte[] Footer = Encoding.ASCII.GetBytes("RAWDONE1");

    private static int Main(string[] args)
    {
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
        CultureInfo.CurrentUICulture = CultureInfo.InvariantCulture;
        if (args.Length == 1 && args[0] == "--version")
        {
            Console.WriteLine("RadiantThermoRawReader 1.0.0; binary protocol 2");
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
        using var raw = RawFileReaderFactory.ReadFile(path);
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

        using var buffer = new BufferedStream(Console.OpenStandardOutput(), 1024 * 1024);
        using var writer = new BinaryWriter(buffer, Encoding.UTF8, leaveOpen: true);
        writer.Write(Magic);
        writer.Write(ProtocolVersion);
        writer.Write(UnknownSpectrumCount);
        uint sequentialNumber = 0;
        for (int nativeNumber = raw.RunHeaderEx.FirstSpectrum;
             nativeNumber <= raw.RunHeaderEx.LastSpectrum;
             ++nativeNumber)
        {
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
                double? width = reaction.IsolationWidth >= 0 ? reaction.IsolationWidth : trailer.AsDouble("MS2 Isolation Width:");
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

            using var frame = new MemoryStream(checked(40 + masses.Length * 8));
            using var frameWriter = new BinaryWriter(frame, Encoding.UTF8, leaveOpen: true);
            frameWriter.Write(checked(36u + (uint)masses.Length * 8u));
            frameWriter.Write(checked((int)++sequentialNumber));
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
            frameWriter.Flush();
            writer.Write(frame.GetBuffer(), 0, checked((int)frame.Length));
        }
        if (sequentialNumber == 0)
            throw new IOException("The RAW file contains no MS1 or MS2 spectra.");

        writer.Write(Footer);
        writer.Write(sequentialNumber);
        writer.Flush();
    }

    private static (double[] Masses, double[] Intensities) ReadCentroids(
        IRawData raw, IScanEvent scanEvent, int scanNumber)
    {
        // SpectrumWriter.ReadMZData is private protected in the parser library.
        // This is its native-centroid selection and array ordering, using the
        // identical RawFileReader API. No centroid algorithm is implemented here.
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
