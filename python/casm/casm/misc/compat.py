"""Functions that are Python2/3 compatible"""
import io
import six

# Notes on Python2/3 compatible io for casm:

# For pandas csv files:
# with open('test.csv', compat.pandas_wmode()) as f:
#     f.write('# ')  # will make this optional in a future version
#     df.to_csv(f, sep=compat.str(' '), index=False)
#
# with open('test.csv', compat.pandas_rmode()) as f:
#     if compat.peek(f) == '#':
#         f.read(1)
#     df = pandas.read_csv(f, sep=compat.str(' +'), engine='python')
#     print(df)

# For json files: 
# with open('test.json', 'wb') as f:
#     f.write(six.u(json.dumps(data, indent=2)).encode('utf-8'))
#
# with open('test.json', 'rb') as f:
#     data = json.loads(f.read().decode('utf-8'))
#     print(data)

# For pickle files:
# with open('test.pkl', 'wb') as f:
#     pickle.dump(data, f, protocol=2)
#
# with open('test.pkl', 'rb') as f:
#     data = pickle.load(f)
#     print(data)

def str(val):
    """Convert to text to Python2/3 native string literal type"""
    if six.PY2:
        return val.encode('utf-8')
    else:
        return val

def peek(f, size=1):
    pos = f.tell()
    data = f.read(size)
    f.seek(pos)
    return data

def pandas_wmode():
    if six.PY2:
        return 'wb'
    else:
        return 'w'

def pandas_rmode():
    return 'r'

def native_io():
    """Return Python2/3 IO stream (BytesIO or StringIO)"""
    if six.PY2:
        return io.BytesIO()
    else:
        return io.StringIO()

