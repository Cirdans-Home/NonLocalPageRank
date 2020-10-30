The data collected in 

london_tube_nodes
london_tube_edges
london_tube_layers

describes a 13-layer multiplex network whose nodes are London train stations, edges are rail connections between the stations and each layer indicates which type of rail system each edge belongs to: overground (layer 12), DLR (layer 13) or one of the 11 underground rail lines (layers 1-11). The overall number of stations (nodes) is 369.

The dataset uses the three-layer multiplex network from

[1] 
M. De Domenico, A. Solé-Ribalta, S. Gómez, and A. Arenas. Navigability of interconnected networks under random failures.Proceedings of the National Academy of Sciences, 111(23):8351–8356, 2014.

as a baseline.

The data collected in 

london_tube_usage

Annual usage of tube stations in million of passengers for 10 years

reports annual passenger usage of underground rail stations in million of passengers over the span of 10 years. As it only considers underground rail stations (the first 11 layers of the multiplex) the dataset consists of only 271 stations (nodes). 
