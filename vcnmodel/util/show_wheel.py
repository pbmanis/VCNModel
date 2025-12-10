import pprint
from zipfile import ZipFile

path = '/Users/pbmanis/Desktop/Python/cnmodel/dist/cnmodel-0.58.2-py3-none-any.whl'
names = ZipFile(path).namelist()
pprint.pprint(names)
