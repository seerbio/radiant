import numpy as np
import pandas as pd
from matplotlib import pyplot as p
import struct

NEUTRON_ISODIFF = 1.0031
HYDROGEN = 1.0078
PROTON = 1.0072
H2O = 18.010565
OXYGEN = 15.99491463
CARBON = 12.00
NITROGEN = 14.003074
AMMONIA = NITROGEN + (3 * HYDROGEN)

AA = dict(A=71.037114,
          C=103.00919,
          D=115.02694,
          E=129.04260,
          F=147.06842,
          G=57.021464,
          H=137.05891,
          I=113.08406,
          K=128.09496,
          L=113.08406,
          M=131.04048,
          N=114.04293,
          P=97.05858,
          Q=128.05858,
          R=156.10111,
          S=87.032029,
          T=101.04768,
          V=99.068414,
          W=186.07931,
          Y=163.06333)


AA_index = dict(D=1, E=2, P=3, M=4, C=5,
                N=6, Q=7,
                S=8, T=9,
                G=10, A=11,
                I=12, L=13, V=14,
                F=15, W=16, Y=17,
                H=18, K=19, R=20)


def peptide_sequence_to_array(seq: str, seq_max_len: int) -> np.array:

    if not len(seq) <= seq_max_len:
        print(seq, len(seq), seq_max_len)

    assert len(seq) <= seq_max_len

    arr = np.zeros(seq_max_len)

    for i, aa, in enumerate(seq):
        arr[i] = AA_index[aa]

    return arr.astype(int)


def to_one_hot(arr: np.array, seq_max_len):

    assert len(arr) <= seq_max_len

    number_of_amino_acids_plus_null = 21
    mat = np.zeros([number_of_amino_acids_plus_null, seq_max_len])

    for i, aa_pos in enumerate(arr):
        mat[aa_pos, i] = 1

    return mat


def deserialize_doubles(serialized_doubles: bytes) -> np.array:
    d_amount = ('d' * int(len(serialized_doubles) / 8))
    return np.array(struct.unpack(d_amount, serialized_doubles))


def deserialize_floats(serialized_floats: bytes) -> np.array:
    f_amount = ('f' * int(len(serialized_floats) / 4))
    return np.array(struct.unpack(f_amount, serialized_floats))


def deserialize_ints(serialized_ints: bytes) -> np.array:
    i_amount = ('i' * int(len(serialized_ints) / 4))
    return np.array(struct.unpack(i_amount, serialized_ints))


def serialize_doubles(doubles_arr: np.array) -> bytes:
    d_amount = ('d' * len(doubles_arr))
    return struct.pack(d_amount, *doubles_arr.tolist())


def serialize_floats(floats_arr: np.array) -> bytes:
    f_amount = ('f' * len(floats_arr))
    return struct.pack(f_amount, *floats_arr.tolist())


def get_max_prediction_length_by_charge(charge: int) -> int:

    if charge == 2:
        return 40
    elif charge > 2:
        return 50

    return 22


def amino_acid_masses(fixed_modifications_dict: dict):

    aa_acid_masses = AA.copy()
    for residue_letter, modification_mass in fixed_modifications_dict.items():
        aa_acid_masses[residue_letter] += modification_mass

    return aa_acid_masses


def calucluate_precursor_mass(seq: str, charge: int, amino_acid_mass: dict) -> float:

    mass = H2O
    for i, aa in enumerate(seq):
        mass += amino_acid_mass[aa]

    return ((charge * PROTON) + mass) / charge


def build_mz_values_list(seq: str,
                         charge: int,
                         start_mass: float,
                         max_length,
                         amino_acid_mass: dict) -> list:

    null_value = -1

    arr = np.zeros(max_length)
    for i, aa in enumerate(seq):

        if i >= len(arr):
            break

        arr[i] = amino_acid_mass[aa]

    arr = np.cumsum(arr) + start_mass
    arr /= charge
    arr[len(seq):] = null_value

    return arr.tolist()


def residue_lengths(charge):
    if charge == 2:
        l1 = 40
        l2 = 40
        a = 40
    elif charge == 1:
        l1 = 22
        l2 = -1
        a = 22
    else:
        l1 = 50
        l2 = 50
        a = 50

    return dict(l1=l1, l2=l2, a=a)


def build_mz_values_new(seq, charge, fixed_modifications_dict, length, length2, a_len):

    max_peptide_length = get_max_prediction_length_by_charge(charge)
    seq_reversed = seq[::-1]
    aa_masses: dict = amino_acid_masses(fixed_modifications_dict)
    precursor_mass = calucluate_precursor_mass(seq, charge, aa_masses)

    if len(seq) > max_peptide_length:
        raise ValueError

    ions = []

    y_series = build_mz_values_list(seq_reversed, 1, HYDROGEN + H2O, length, aa_masses)
    a_series = build_mz_values_list(seq, 1, HYDROGEN - (CARBON + OXYGEN), a_len, aa_masses)
    b_series = build_mz_values_list(seq, 1, HYDROGEN, length, aa_masses)

    if charge < 3:
        ions.append(precursor_mass)
        ions.append(precursor_mass - H2O)
        ions.append(precursor_mass - AMMONIA)

    if len(seq) <= length:
        b_series[len(seq) - 1] += H2O

    ions += y_series

    if charge > 1:
        y_charge_2_series = build_mz_values_list(seq_reversed, 2, (2 * HYDROGEN) + H2O, length2, aa_masses)
        ions += y_charge_2_series

    y_nh3_series = np.array(y_series) - AMMONIA
    y_nh3_series[y_nh3_series < 0] = -1.0
    ions += list(y_nh3_series)

    y_h2o_series = np.array(y_series) - H2O
    y_h2o_series[y_h2o_series < 0] = -1.0
    ions += list(y_h2o_series)

    ions += b_series
    if charge > 1:
        b_charge_2_series = build_mz_values_list(seq, 2, (2 * HYDROGEN), length2, aa_masses)
        if len(seq) <= length2:
            b_charge_2_series[len(seq) - 1] += H2O / 2
        ions += b_charge_2_series

    b_nh3_series = np.array(b_series) - AMMONIA
    b_nh3_series[b_nh3_series < 0] = -1.0
    ions += list(b_nh3_series)

    b_h2o_series = np.array(b_series) - H2O
    b_h2o_series[b_h2o_series < 0] = -1.0
    ions += list(b_h2o_series)

    ions += a_series

    return ions


def build_frags(charge):

    res_len = residue_lengths(charge)

    length = int(res_len['l1'])
    length2 = int(res_len['l2'])
    a_len = int(res_len['a'])

    cols = []
    for x in range(1, length + 1):

        cols.append('y%s' % x)

    for x in range(1, length + 1):
        cols.append('b%s' % x)

    for x in range(1, a_len + 1):
        cols.append('a%s' % x)

    if charge > 1:
        for x in range(1, length2 + 1):
            cols.append('y%s^2' % x)

        for x in range(1, length2 + 1):
            cols.append('b%s^2' % x)

    for x in range(1, length + 1):
        cols.append('y%s-NH3' % x)

    for x in range(1, length + 1):
        cols.append('b%s-NH3' % x)

    for x in range(1, length + 1):
        cols.append('y%s-H2O' % x)

    for x in range(1, length + 1):
        cols.append('b%s-H2O' % x)

    return cols


def plot_msms_raw(mz, intensity, title, annotate_fraction: float = 2):

    intensity /= intensity.max()

    for m, i in zip(mz,intensity):
        p.plot([m, m], [0, i], 'b')

        if (annotate_fraction < 1) and (i > annotate_fraction):
            p.text(m, i + 0.01, round(m, 2), rotation=90)

    p.xlabel("m/z")
    p.ylabel("Rel Intensity")
    p.axhline(y=0, color='b', linestyle='-')
    p.title(title)

    p.show()


def mirror_plot_msms_raw(mz1, intensity1, mz2, intensity2, title: str = "", labels=None, annotate_fraction: float = 2):

    intensity1 /= intensity1.max()

    if labels is None:
        labels = list()
    p.axhline(y=0, color='b', linestyle='-')
    p.title(title)
    for m, i in zip(mz1, intensity1):
        p.plot([m, m], [0, i], 'b')

    for m, i in zip(mz2, intensity2):
        p.plot([m, m], [0, -i], 'b')

    if len(labels) > 0:
        for m, i, l in zip(mz2, intensity2, labels):
            if (annotate_fraction < 1) & (-i < -annotate_fraction):
                p.text(m, -i - 0.1, l, rotation=90)

    p.xlabel("m/z")
    p.ylabel("Rel Intensity")

    p.show()


def _generate_ions_by_type(max_len: int, ion_type: str, modifier: str = ""):
    return [f'{ion_type}{x}{modifier}' for x in range(1, max_len + 1)]


def generate_all_ions(charge, seq_max_length_charge):

    ions_to_generate = [('y', ""), ('y', "^2"), ('y', "-NH3"), ('y', "-H2O"), ('b', ""), ('b', "^2"), ('b', "-NH3"),
                        ('b', "-H2O"), ('a', "")]

    generated_ions = []
    if charge == 1:
        del ions_to_generate[5]
        del ions_to_generate[1]
        generated_ions += ['p', 'p-H2O', 'p-NH3']

    if charge == 2:
        generated_ions += ['p^2', 'p-H2O^2', 'p-NH3^2']

    for itg in ions_to_generate:

        ions = _generate_ions_by_type(seq_max_length_charge, itg[0], itg[1])
        generated_ions += ions

    return generated_ions


def add_missing_columns(df: pd.DataFrame, ions_in_model):

    df_cols = df.columns
    for ion in ions_in_model:
        if ion not in df_cols:
            df[ion] = 0

    return df


def get_max_seq_length_by_charge(charge):

    if charge == 2:
        return 40
    elif charge > 2:
        return 50
    else:
        return 22


def calc_normalized_collision_energy(collision_energy_ev: int, charge, precursor_mz_thomsons):

    thermo_mz_normalizer: float = 500.0
    charge_factor: float = 1.0

    if charge == 2:
        charge_factor = 0.9
    elif charge == 3:
        charge_factor = 0.85
    elif charge == 4:
        charge_factor = 0.80

    return collision_energy_ev / ((precursor_mz_thomsons / thermo_mz_normalizer) * charge_factor)
