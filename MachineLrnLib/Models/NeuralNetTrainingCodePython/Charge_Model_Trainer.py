import pandas as pd
import numpy as np
import random

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

rando = 666

####*IMPORANT*: Have to do this line *before* importing tensorflow
import os
os.environ['PYTHONHASHSEED']=str(rando)

from keras.models import Sequential
from keras.layers import Dense
from keras.layers import BatchNormalization
from keras.layers import Flatten
from keras.layers import Activation
from keras.callbacks import ModelCheckpoint
from keras.layers import Conv1D
from keras.layers import Dropout
from keras.layers import LSTM
from keras.layers import Bidirectional
from keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
# from keras.utils import to_categorical
from keras import utils as np_utils
import utils as ut
# import tensorflow as tf
import multiprocessing as mp
from multiprocessing import freeze_support
from tqdm import tqdm
import functools
import time
import convert_model as cm
# tf.random.set_seed(rando)
random.seed(rando)
np.random.seed(rando)
from sklearn.preprocessing import Normalizer


charge = 4
seq_max_length_charge = ut.get_max_seq_length_by_charge(charge)


def baseline_model_rnn(X, Y):

    model = Sequential()
    num_classes = Y.shape[1]

    nodes = -1
    if charge == 1:
        nodes = 48
    if charge == 2:
        nodes = 128
    if charge >= 3:
        nodes = 192

    model.add(Bidirectional(LSTM(nodes, return_sequences=True, activation='tanh'), input_shape=(X.shape[1:])))
    model.add(Bidirectional(LSTM(nodes, activation='tanh')))
    # model.add(Dropout(0.1))
    model.add(Dense(num_classes, activation='linear', use_bias=True))
    model.compile(loss='mean_squared_error',
                  optimizer='adam',
                  metrics=['mean_squared_error', 'cosine_proximity'])

    return model


def build_x_data(seqs: list, seq_max_length: int) -> list:

    pool = mp.Pool(mp.cpu_count())
    results = pool.map(functools.partial(ut.peptide_sequence_to_array, seq_max_len=seq_max_length), seqs)

    x_one_hots = []
    for r in results:
        x_one_hots.append(ut.to_one_hot(r, seq_max_length))

    pool.close()

    return x_one_hots


if __name__ == '__main__':

    freeze_support()

    filepath = r"C:\Repositories\PyCharm\FragPredData\all_combined_nce_%s.csv" % charge
    print(filepath)
    df = pd.read_csv(filepath)
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
    print(time.time() - t)

    for x, i in zip(x_vals, nce):
        x *= i

    loss = 'loss'
    callbacks = [
        EarlyStopping(patience=10, verbose=1, monitor=loss),
        ReduceLROnPlateau(factor=0.1, patience=5, min_lr=0.00001, verbose=1, monitor=loss),
        ModelCheckpoint('checkpoints2.h5', verbose=False, save_best_only=True, save_weights_only=True, monitor=loss)
    ]

    if charge == 1:
        batch_size = 20
    elif charge == 2:
        batch_size = 400
    elif charge == 3:
        batch_size = 250
    else:
        batch_size = 100

    model = baseline_model_rnn(x_vals, y_vals)
    model.fit(x_vals, y_vals, validation_data=(x_vals, y_vals), epochs=100, batch_size=batch_size, callbacks=[callbacks])

    model_filepath = "Models/rnn_linear_charge_w_precursors_nce_%s.hdf5" % charge
    model.save(model_filepath)

    y_cols = ','.join(filtered_cols_y)
    print(y_cols)
