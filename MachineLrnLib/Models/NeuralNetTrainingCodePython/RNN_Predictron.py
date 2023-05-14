import numpy as np
import os
import pandas as pd
import utils as ut

from keras.models import load_model

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

np.set_printoptions(linewidth=4000)

# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


class PredictronRNN:

    def __init__(self,
                 charge: int,
                 fixed_modifications_dict: dict,
                 intensity_threshold=0.01,
                 serialize_arrays: bool = True):

        self.charge: int = charge
        self.neural_network_model: str = self._get_model_filepath()
        self.model = load_model(self.neural_network_model)
        self.max_length = ut.get_max_prediction_length_by_charge(charge)
        self.intensity_threshold = intensity_threshold
        self.ion_cols = ut.generate_all_ions(self.charge, self.max_length)
        self.fixed_modifications_dict = fixed_modifications_dict
        self.serialize_arrays = serialize_arrays

    def predict_msms(self, list_of_sequences: list, list_of_norm_coll_energy) -> pd.DataFrame:

        x_list: list = []
        mz_list: list = []
        mz_new: list = []
        seq_list: list = []

        res_lens = ut.residue_lengths(self.charge)

        for seq, nce in zip(list_of_sequences, list_of_norm_coll_energy):

            if (len(seq) > self.max_length) | (not self._sequence_is_valid(seq)):
                continue

            seq_list.append(seq)

            mat = ut.peptide_sequence_to_array(seq, self.max_length)
            x_list.append(ut.to_one_hot(mat, self.max_length) * nce)

            mz = ut.build_mz_values_new(seq, self.charge, self.fixed_modifications_dict,
                                         res_lens['l1'], res_lens['l2'], res_lens['a'])

            mz_list.append(mz)

        intensities = self.model.predict(np.array(x_list))
        intensities = np.round((intensities / intensities.max()), 2)

        assert len(intensities) == len(mz_list)
        assert len(intensities) == len(seq_list)
        assert len(intensities[0]) == len(self.ion_cols)
        assert len(mz_list[0]) == len(self.ion_cols)

        for mz, intensity in zip(mz_list, intensities):
            df_seq = pd.DataFrame()
            df_seq['mz'] = mz
            df_seq['intensity'] = intensity
            df_seq['ions'] = self.ion_cols

            df_seq = df_seq.loc[(df_seq.mz > 0) & (df_seq.intensity >= self.intensity_threshold), :]
            mz_new.append(df_seq)

        df: pd.DataFrame = pd.DataFrame()
        df['sequence'] = seq_list
        df['prediction'] = mz_new
        return df

    def _get_model_filepath(self):
        assert (self.charge > 0) & (self.charge < 5)
        # return "Models/rnn_tanh_charge_%s.hdf5" % self.charge
        return "Models/rnn_linear_charge_w_precursors_nce_%s.hdf5" % self.charge

    def _sequence_is_valid(self, seq: str):

        invalid_amino_acid_letters = ['B', 'J', 'O', 'U', 'X', 'Z']
        seq = seq.upper()

        for invalid_amino_acid_letter in invalid_amino_acid_letters:

            if invalid_amino_acid_letter in seq:
                return False

        return True


if False:
    charge = 2

    # seq = 'DTLMISR'
    # seq = 'EFHLHLR'
    # seq = 'QFHLHLR'
    seq = 'ALDEDGK'
    # seq = 'ALDGEDK'
    seq = 'GKFLLSLAYSPDGK'
    seq = 'SAMVEDLSLR'
    # seq = 'QHFLHLR'
    seq2 = 'AMAVEDLLSR'

    mods = dict(C=57.021464)

    pred = PredictronRNN(charge, mods)

    import time

    t = time.time()
    preds = pred.predict_msms([seq, seq2], [0.30, 0.30])
    pr = preds.iloc[0]
    pr2 = preds.iloc[1]
    print(time.time() - t)
    print(pr2.prediction)

    points = [
        [120.081, 43709.5],
        [129.09, 7363.5],
        [129.102, 182194],
        [129.112, 4607.4],
        [130.065, 49165.5],
        [130.087, 55051.5],
        [131.069, 6989.18],
        [132.081, 20052.3],
        [136.076, 43849],
        [136.589, 18084.2],
        [143.082, 45527.9],
        [143.118, 1.04172e+06],
        [144.08, 7278.2],
        [144.121, 9395.83],
        [146.06, 33964.4],
        [147.113, 121665],
        [148.116, 5352.33],
        [150.055, 7398.26],
        [158.092, 6588.98],
        [159.075, 13113.4],
        [159.092, 351032],
        [160.078, 7346.48],
        [160.095, 8479.34],
        [167.082, 21130.1],
        [169.077, 9479.23],
        [170.06, 29084.8],
        [171.078, 6836.71],
        [171.113, 117839],
        [171.131, 9013.78],
        [173.073, 7635.47],
        [173.129, 7314.57],
        [174.055, 19815.4],
        [175.119, 46002.8],
        [182.363, 5114],
        [183.113, 8224.42],
        [185.128, 6948.3],
        [187.087, 68244.8],
        [189.067, 5723.53],
        [191.085, 95070.5],
        [194.094, 7338.97],
        [195.076, 30934.2],
        [200.14, 8299.15],
        [201.067, 6469.33],
        [201.123, 13628.7],
        [202.05, 14663.1],
        [203.16, 5594.69],
        [204.135, 7527.7],
        [209.092, 17041.4],
        [215.14, 8209.42],
        [217.138, 7219.27],
        [219.08, 141443],
        [225.101, 7483.17],
        [226.119, 47775.8],
        [227.285, 36174.6],
        [228.134, 35794],
        [229.119, 7870.96],
        [240.134, 14351.3],
        [243.147, 25976],
        [244.13, 61108.6],
        [252.114, 15896.4],
        [258.146, 15409],
        [261.121, 7660.62],
        [261.156, 215710],
        [270.125, 15782.2],
        [271.027, 6899.76],
        [272.137, 9292.79],
        [272.172, 43143.9],
        [287.145, 7366.66],
        [294.145, 5382.08],
        [295.177, 6832.17],
        [317.183, 11890.8],
        [320.102, 6765.54],
        [332.166, 45108.5],
        [335.153, 11150.9],
        [336.138, 6422.42],
        [337.123, 10982.8],
        [339.204, 20511.3],
        [344.197, 8964.21],
        [353.16, 7487.56],
        [357.217, 26643.2],
        [359.173, 7728.94],
        [363.146, 8542.74],
        [365.186, 9743.92],
        [372.63, 23020.9],
        [372.706, 10377.1],
        [374.24, 146397],
        [377.491, 13264.9],
        [377.821, 14443.6],
        [378.152, 8427.81],
        [381.155, 87214.5],
        [381.267, 10879.5],
        [401.186, 10003.5],
        [405.21, 18935.8],
        [407.652, 5751.86],
        [407.82, 13357.5],
        [408.158, 6406.22],
        [415.661, 13369.3],
        [446.208, 14406.2],
        [447.193, 7735.37],
        [453.15, 7465.91],
        [461.168, 21051.5],
        [461.665, 7547.53],
        [463.247, 9050.56],
        [480.25, 6788.86],
        [505.281, 62780.6],
        [518.297, 11847.7],
        [542.695, 9776.5],
        [544.248, 11355.3],
        [574.3, 9543.75],
        [575.285, 20504.2],
        [592.312, 440062],
        [593.315, 18044.3],
        [602.296, 8882.07],
        [610.322, 6578.73],
        [637.142, 6319.81]
    ]

    rr = np.array(points)

    # ut.mirror_plot_msms_raw(rr[:,0], rr[:,1],
    #                         pr.prediction.mz, pr.prediction.intensity,
    #                         seq, pr.prediction.ions, 0.02)

    ut.mirror_plot_msms_raw(pr2.prediction.mz, pr2.prediction.intensity,
                            pr.prediction.mz, pr.prediction.intensity,
                            seq, pr.prediction.ions, 0.02)

    # filepath = r"C:\Repositories\PyCharm\FragPredData\all_combined_nce_%s.csv" % charge
    # df = pd.read_csv(filepath)
    # df = df.loc[df.Sequence == seq]
    # print(df)
    # df = df.loc[:, pr.prediction.ions.values]
    #
    # pred = df.iloc[0].values
    # pred_int = pred / pred.max()
    #
    # print(pred)

    # ut.mirror_plot_msms_raw(pr.prediction.mz, pr.prediction.intensity, pr.prediction.mz, pred_int)
