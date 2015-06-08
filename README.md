# SACLA2015A8038
Data analysis repository for SACLA experiment between June 8-12, 2015 (proposal number: 2015A8038)

--------------------------------------------------------------------------------

DataConvert
-----------

- Login to one of the analysis nodes at SACLA:

  `ssh fperakis@xhpcfep`
  `qsub -I`

- Run `DataConvert3` to extract the data to HDF5 format:

  `DataConvert3 -dir ~/data/test -o 259782.h5 -each -mode 30 -r 259782 -starttag 447317100 -endtag 447352620`

- Run `DataCompress3.py` using Python 2.7 to produce averages for each scan position:

  `python2.7 dataCompress3.py -f 259782.h5 -c run259782_out.csv -a`

