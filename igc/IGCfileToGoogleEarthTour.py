import os
import numpy as np
import pandas as pd
import simplekml

# Import CSV file from current directory
path = '2023-02-22-XNA-9AD179D6D6F04DDD95474658245E6D4F-01'+'.csv'
df_B = pd.read_csv(os.path.join(os.getcwd(),path),header=0,index_col=0)

# Convert dataframe to list of tuple
tuple_B =  [tuple(x) for x in df_B[['longitude','latitude','press_alt']].values]

# Convert longitude/latitude to the distance from the point of beginning
dt = 1
lat = np.radians(df_B.latitude.values-df_B.latitude.values[0])*6378137
lng = np.radians(df_B.longitude.values-df_B.longitude.values[0])*6378137
alt = df_B.press_alt.values

# Calculate the angles of yaw (heading), pitch (tilt) and roll from the coordinate
df_B['heading'] = 90
df_B['tilt'] = 90
df_B['roll'] = 0
for n in np.arange(len(df_B)):
  if dt-1 < n and n < len(df_B)-dt:
    a = lat[n+dt]-lat[n-dt]
    b = lng[n+dt]-lng[n-dt]
    c = lng[n-dt]*lat[n+dt]-lng[n+dt]*lat[n-dt]
    d = np.sqrt(a*a+b*b)
    e = np.abs(-a*lng[n]+b*lat[n]+c)/d
    # Calculate turn radius and rate
    R = (d*d+4*e*e)/(8*e)
    Omega = np.arcsin((d/2)/R)/dt
    # Calculate yaw, pitch and roll
    df_B.loc[n,'heading'] = -np.rad2deg(np.arctan2(a,b))-(-90)
    df_B.loc[n,'tilt'] = 90+np.rad2deg(np.arctan2(alt[n+dt]-alt[n-dt],d))
    df_B.loc[n,'roll'] = np.rad2deg(np.arctan2(R*Omega*Omega,9.81))
    # Modification
    df_B.loc[n,'heading'] = df_B.heading[n]-360 if df_B.heading[n] > 180 else df_B.heading[n]
    df_B.loc[n,'roll'] = -df_B.roll[n] if -a*lng[n]+b*lat[n]+c < 0 else df_B.roll[n]

# Insert values to beginning and end of dataframe
df_B.loc[0:dt,'heading'] = df_B.heading[dt]
df_B.loc[0:dt,'tilt'] = df_B.tilt[dt]
df_B.loc[len(df_B)-dt:,'heading'] = df_B.heading[len(df_B)-dt-1]
df_B.loc[len(df_B)-dt:,'tilt'] = df_B.tilt[len(df_B)-dt-1]

## Create an instance of Kml
kml = simplekml.Kml(open=1)

# Create a linestring
linestring = kml.newlinestring(name="A Sloped Line")
linestring.coords = tuple_B
linestring.altitudemode = simplekml.AltitudeMode.relativetoground
linestring.extrude = 0
linestring.style.linestyle.width = 3
linestring.style.linestyle.color = simplekml.Color.red

# Create a tour and attach a playlist to it
tour = kml.newgxtour(name="Take-Off!")
playlist = tour.newgxplaylist()

n=0
for _ in df_B.iterrows():
    # Attach a gx:FlyTo to the playlist
    flyto = playlist.newgxflyto()
    flyto.camera.longitude = df_B.longitude[n]
    flyto.camera.latitude = df_B.latitude[n]
    flyto.camera.altitude = df_B.press_alt[n]+1.
    flyto.camera.heading = df_B.heading[n]
    flyto.camera.tilt = df_B.tilt[n]
    flyto.camera.roll = df_B.roll[n]
    flyto.gxduration = 1.
    flyto.gxflytomode = 'smooth'
    n+=1

# Attach a gx:Wait to the playlist to give the gx:AnimatedUpdate time to finish
wait = playlist.newgxwait(gxduration=3)

# Save to file
kml.save(os.path.join(os.getcwd(),path[:-4]+".kml"))