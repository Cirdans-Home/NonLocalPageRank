%% london_geomap.m
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute

shpfile =  "London/greater-london-latest-free.shp/gis_osm_roads_free_1";
S = shaperead(shpfile);
mapshow(S);