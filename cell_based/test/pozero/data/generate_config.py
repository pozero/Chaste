v0 = [
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
]

p0 = [
    3.0,
    3.1,
    3.2,
    3.3,
    3.4,
    3.5,
    3.6,
    3.7,
    3.8,
    3.9,
    4.0,
]

name_format = '{condition}_{p:d}_{v:02d}'
file_format = '''random_propelling_parameter: {v}
uniform_taret_area: {b:d}
shape_index: {p}
'''

for con in ['FinalNonuniform', 'FinalUniform']:
    for p in p0:
        for v in v0:
            file_name = name_format.format(condition = con, p = int(p * 10), v = v)
            f = open(file_name, 'w')
            f.write(file_format.format(v = v, b = (con == 'FinalUniform'), p = p))
            f.close()