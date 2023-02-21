import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import IGC file from current directory
path = '2023-02-22-XNA-9AD179D6D6F04DDD95474658245E6D4F-01.igc'
with open(path) as f:
    igcfile = f.readlines()

# Extract B-record from IGC file
date = datetime.date(int(path.split('-')[0]),int(path.split('-')[1]),int(path.split('-')[2]))
igcfile = [line.rstrip("\n") for line in igcfile]
igc_B = [[igcfile.index(lines),lines] for lines in igcfile if 'B' in lines[0] if 'A' in lines]
list_B = [[lines[0],lines[1][0],lines[1][1:7],lines[1][7:14],lines[1][14],lines[1][15:23],lines[1][23],lines[1][24],lines[1][25:30],lines[1][30:]] for lines in igc_B]

# Convert text of B-record to values
idx=0
for _ in list_B:
  list_B[idx][2] = datetime.datetime(date.year,date.month,date.day,int(list_B[idx][2][:2]),int(list_B[idx][2][2:4]),int(list_B[idx][2][4:]))
  list_B[idx][3] = float(list_B[idx][3])/100000.
  list_B[idx][5] = float(list_B[idx][5])/100000.
  list_B[idx][8] = float(list_B[idx][8])
  list_B[idx][9] = float(list_B[idx][9])
  idx+=1

# Create dataframe from list of records and export it as CSV file
df_B = pd.DataFrame(list_B,columns=['igc_index', 'record_type', 'UTC_time','latitude','NS','longitude','EW','fix_validity','press_alt','GNSS_alt'])
df_B.to_csv(path[:-4]+".csv")

df_B0 = df_B.copy()

# Low pass filtering by using FFT
for column in ['latitude','longitude','press_alt','GNSS_alt']:
  f = df_B[column].values
  F = np.fft.rfft(f)
  freq = np.fft.rfftfreq(len(f))
  F[(1/(freq+1e-10) < 20)] = 0
  df_B[column] = np.fft.irfft(F, len(f))

# Convert dataframe to list in same format as IGC
list_B_out = [[lines[0],lines[1]+lines[2].strftime('%H%M%S')+format(int(lines[3]*100000),'07')+lines[4]+format(int(lines[5]*100000),'08')+lines[6]+lines[7]+format(int(lines[8]),'05')+format(int(lines[9]),'05')] for lines in df_B.values.tolist()]
igcfile_out = igcfile
idx_df=0
for idx in np.arange(len(igcfile_out)):
  if idx <= df_B.igc_index.max() and df_B.values.tolist()[idx_df][0] == idx:
    line_B = df_B.values.tolist()[idx_df]
    igcfile_out[idx] = line_B[1]+line_B[2].strftime('%H%M%S')+format(int(line_B[3]*100000),'07')+line_B[4]+format(int(line_B[5]*100000),'08')+line_B[6]+line_B[7]+format(int(line_B[8]),'05')+format(int(line_B[9]),'05')
    idx_df+=1

# Save created list as new igc file
path_w = path[:-4]+"FFT"+path[-4:]
with open(path_w, mode='w') as f:
    f.write('\n'.join(igcfile_out))

# Plot 3D flight path
fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(projection='3d')

x = np.radians(df_B.longitude-df_B.longitude[0])*6378.137*1000
y = np.radians(df_B.latitude-df_B.latitude[0])*6378.137*1000
z = df_B.press_alt-df_B.press_alt.min()
range = max(abs(x.min()),abs(x.max()),abs(y.min()),abs(y.max()),abs(z.min()),abs(z.max()))

ax.view_init(elev=30, azim=-60)
ax.set_xlim3d(-range,range)
ax.set_ylim3d(-range,range)
ax.set_zlim3d(0,range)
ax.set_box_aspect((2,2,1))
ax.plot(x,y,z,c="red")
plt.show()
