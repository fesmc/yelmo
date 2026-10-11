# ytrc

| id | variable          | dimensions       | units       | long_name                                            |
|----|-------------------|------------------|-------------|------------------------------------------------------|
|  1 | t_dep             | xc, yc, zeta     | yr          | Deposition time, authoritative backend               |
|  2 | t_dep_euler       | xc, yc, zeta     | yr          | Deposition time, Eulerian backend                    |
|  3 | t_dep_euler_rel   | xc, yc, zeta     | yr          | Deposition time minus t_ref_euler, Eulerian backend  |
|  4 | t_dep_trc         | xc, yc, zeta     | yr          | Deposition time, Lagrangian particle backend         |
|  5 | t_dep_elsa        | xc, yc, zeta     | yr          | Deposition time, Lagrangian layer (elsa) backend     |
|  6 | depth_iso         | xc, yc, time_iso | m           | Depth of isochronal layers, authoritative backend    |
|  7 | trc_count         | xc, yc, depth_norm | 1         | Tracer particle count per normalized-depth band      |
|  8 | trc_depth_iso     | xc, yc, time_iso | m           | Isochrone depth from the particle cloud (tracer)     |
