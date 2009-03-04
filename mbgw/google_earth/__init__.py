# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

str_1 = r"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.2">
  <Document>
    <name>PfPR</name>
    <description>Age-adjusted PfPR, 2-10 year old age group</description>
    <Style id="highlightPlacemark">
      # <IconStyle>
      #   <Icon>
      #     <href>http://maps.google.com/mapfiles/kml/paddle/red-stars.png</href>
      #   </Icon>
      # </IconStyle>
      <BalloonStyle>
        <text>
          <![CDATA[
            Latitude: $[Lat]
            Longitude: $[Lon]
            Number examined: $[Examined]
            Number positive: $[Pf]
            Parasite rate: $[PfPR]
            Age-corrected parasite rate $[PfPR210]
          ]]>
        </text>
      </BalloonStyle>      
    </Style>
    <Style id="redPlacemark">
      <IconStyle>
        <Icon>
          # <href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href>
          <href>red.png</href>
        </Icon>
      </IconStyle>
    </Style>
    <Style id="orangePlacemark">
      <IconStyle>
        <Icon>
          # <href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href>
          <href>arng.png</href>
        </Icon>
      </IconStyle>
    </Style>
    <Style id="bluePlacemark">
      <IconStyle>
        <Icon>
          # <href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href>
          <href>bloo.png</href>
        </Icon>
      </IconStyle>
    </Style>
    <StyleMap id="exampleStyleMap">
      <Pair>
        <key>normal</key>
        <styleUrl>#normalPlacemark</styleUrl>
      </Pair>
      <Pair>
        <key>highlight</key>
        <styleUrl>#normalPlacemark</styleUrl>
      </Pair>
    </StyleMap>
"""
str_2 = r"""
  </Document>
</kml>
"""

"REC_ID","LAT","LONG","MONTH_STAR","YEAR_START","MONTH_END","YEAR_END","LOW_AGE","UP_AGE","EXAMINED","PF","PFPR","PFPR210"

from csv import reader

def data_to_kml(fnames_in, fname_out):

    f_out = file(fname_out,'w')
    f_out.write(str_1)

    for f_in in fnames_in:

        keys = f_in.next()
        for line in f_in:
        
            line_dict = dict(zip(keys,line))
        
            pfpr = float(line_dict['PFPR210'])
            if pfpr < 5:
                col='blue'
            elif pfpr < 40:
                col='orange'
            else:
                col='red'

            f_out.write("""<Placemark>
              <styleUrl>#%sPlacemark</styleUrl>
    """ % col)            
            f_out.write("\t\t<ExtendedData>\n")
            for key in keys:
                f_out.write("""     <Data name="%s">
            <value>%s</value>
        </Data> 
    """ % (key,line_dict[key]))
            f_out.write("\t\t</ExtendedData>\n")

            f_out.write('\t\t<Point>\n\t\t\t<coordinates>%s,%s,0</coordinates>\n\t\t</Point>\n'    % (line[2],line[1]))
            f_out.write("""    </Placemark>
        """)
    f_out.write(str_2)
    f_out.close()
