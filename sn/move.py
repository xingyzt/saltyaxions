import os

import shutil
masses = [ 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 ]

for m in masses:
    path = f'm{m:04.1f}_g-10.00_0'
    src = 'csv/' + path + '/index.csv'
    dst = f'csv_simple_new/{m:.0f}.csv'
    shutil.copyfile(src, dst)
    
