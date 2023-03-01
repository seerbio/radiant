from keras.models import load_model
import numpy as np
import pandas as pd
import utils as ut
from Charge_Model_Trainer import build_x_data
import time
import multiprocessing as mp
from multiprocessing import freeze_support

def main():

    for charge in range(1, 5):

        filepath_base: str = r"C:/Repositories/PyCharm/FragPredData/mouse_cv_%s.csv" % charge
        df: pd.DataFrame = pd.read_csv(filepath_base)

        seq_max_length_charge = ut.get_max_seq_length_by_charge(charge)

        print(df.shape)
        df = df.loc[df.Sequence_Length <= seq_max_length_charge, :]
        df = df.fillna(0)

        filtered_cols_y = ut.generate_all_ions(charge, seq_max_length_charge)
        df = ut.add_missing_columns(df, filtered_cols_y)

        filtered_cols = ['Sequence', 'NCE']
        filtered_cols += filtered_cols_y

        df = df.loc[:, filtered_cols]
        print(filtered_cols_y)

        df = df.loc[df.iloc[:, 6:].max(axis=1) > 1, :]

        print("Zero NCE Count", df.loc[df.NCE < 1, :].shape)
        df = df.loc[df.NCE > 1, :]

        nce = (df.NCE.values / 100.0)
        print('nce range', nce.min(), nce.max())

        y_vals = df.loc[:, filtered_cols_y]
        print('zeds', np.where(y_vals.max(axis=1).values < 1)[0])
        y_vals = y_vals.values / y_vals.max(axis=1).values.reshape(-1, 1)

        t = time.time()
        x_vals = np.array(build_x_data(df.Sequence, seq_max_length_charge))
        print(x_vals.shape)
        print(filepath_base)
        print(time.time() - t)

        for x, i in zip(x_vals, nce):
            x *= i

        model_filepath: str = "Models/rnn_linear_charge_w_precursors_nce_%s.hdf5" % charge
        model = load_model(model_filepath)

        preds = model.predict(x_vals)
        preds[preds < 0.01] = 0

        from sklearn.metrics import mean_squared_error
        from math import sqrt

        rmse = sqrt(mean_squared_error(y_vals, preds))
        print(charge, 'RMSE', rmse)

        for seq, nc, act, pr in zip(df.Sequence, nce, y_vals, preds):

            res_lens = ut.residue_lengths(charge)
            mz = ut.build_mz_values_new(seq, charge, {},
                                        res_lens['l1'], res_lens['l2'], res_lens['a'])
            print(seq, nc)
            ut.mirror_plot_msms_raw(mz, pr, mz, act)


if __name__ == '__main__':
    freeze_support()
    main()

