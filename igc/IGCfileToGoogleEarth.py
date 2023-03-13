import os
import pandas as pd
import simplekml

# Import CSV file from current directory
path = '2023-02-22-XNA-9AD179D6D6F04DDD95474658245E6D4F-01'+'.csv'
df_B = pd.read_csv(os.path.join(os.getcwd(),path),header=0,index_col=0)

# Convert dataframe to list of tuple
tuple_B =  [tuple(x) for x in df_B[['longitude','latitude','press_alt']].values]

# Create an instance of Kml
kml = simplekml.Kml(open=1)

# Create a linestring
linestring = kml.newlinestring(name="A Sloped Line")
linestring.coords = tuple_B
linestring.altitudemode = simplekml.AltitudeMode.relativetoground
linestring.extrude = 0
linestring.style.linestyle.width = 3
linestring.style.linestyle.color = simplekml.Color.red

# Save to file
kml.save(os.path.join(os.getcwd(),path[:-4]+".kml"))