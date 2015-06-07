# SACLA2015A8038
Data analysis repository for SACLA experiment between June 8-12, 2015 (proposal number: 2015A8038)

--------------------------------------------------------------------------------

DataConvert
-----------

- Run DataConvert3 on one of the analysis nodes of `xhpcfep`:

  `DataConvert3 -dir ~/data/test -o 259782.h5 -each -mode 20 -r 259782 -starttag 447317100 -endtag 447352620`

- Run DataCompress3.py to produce averages for each scan position:

  `./dataCompress3.py -f 259782.h5 -c run259782_out.csv -a`

