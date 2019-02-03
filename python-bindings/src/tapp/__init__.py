import math

from .tapp import *

def tic(raw_data, min_rt = -math.inf, max_rt = math.inf):
    rt = []
    intensity = []
    for i in range(0, len(raw_data.scans)):
        if (raw_data.scans[i].retention_time < min_rt or
            raw_data.scans[i].retention_time > max_rt):
            continue
        sum = 0
        for j in range(0, raw_data.scans[i].num_points):
            sum = sum + raw_data.scans[i].intensity[j]
        intensity = intensity + [sum]
        rt = rt + [raw_data.scans[i].retention_time]

    return (rt, intensity)

def load_example_data():
     raw_data = read_mzxml(
            '/data/qatar/17122018/mzXML/Acute2U_3001.mzXML',
            instrument_type = 'orbitrap',
            resolution_ms1 = 70000,
            resolution_msn = 30000,
            reference_mz = 200,
            polarity = 'pos',
        )
     return raw_data

def example_pipeline():
    print("Loading data")
    raw_data = load_example_data()
    n, m = calculate_dimensions(raw_data, 9, 10, 10)
    print("Estimated memory consumption of the [{0}x{1}] grid: {2:.2f} (MB)".format(n, m, n * m /1024/1024))
    mesh = resample(raw_data, 9, 10, 10)
    return (raw_data, mesh)

RawData.tic = tic
