# UNS2STR, a tool to recover nodes indexing in unstructured quadrilateral datasets.

When a structured mesh is initially generated around an aerofoil shape, indexing information is available on all grid nodes. This information is however lost when the grid is converted into the (unordered) format required by unstructured solvers. Nevertheless, flow-field solution data obtained on such grids can be easily converted back to ordered, recovering known index information from the original mesh file and mapping those data onto the CFD solution.

**UNS2STR** provides this data conversion, taking care of the nodes duplication at the *wake line* necessary to *cut* the unstructured domain, restoring the C-grid topology. The code works on single-block C-shaped quadrilateral grids and can also be applied on multiple CFD solutions is moving grid problems.

From the command line, it is executed by typing:

`UNS2STR str_grid_file uns_sol_file verbosity_flag scale_f`

where the first input argument (after executable name) is the original structured grid file, which is used to provide indexing for all grid nodes. It must be provided as a Tecplot binary load-on-demand file (*.szplt*), containing at x,y,z coordinates of the 2D mesh, lying on the xz plane (i.e. y coordinates are identically zero, but need to be included in the file).

The second argument is the unstructured CFD solution, obtained on that grid once converted to the solver format (unordered). It must also be provided as a Tecplot binary load-on-demand file (*.szplt*), containing coordinates of the unstructured 2D mesh, lying on the xz plane (i.e. y coordinates are identically zero) as well as all the other flor variables stored at grid nodes.

The third input is a flag (0/1) specifying the minimm (0) or maximum (1) vebosity level.

The fourth argument is an optional scaling factor, using during nodes coordinate matching procedure between structured and unstructured datasets. The default value of 1.0 is appropriate for most cases.

**UNS2STR** provides as an output the CFD solution converted into structured (ordered) Tecplot format (default name is "*SU2.szplt*", as well as a topological file (*TOPO.szplt*) which contains the input structured grid along with an additional SU2_NODE-ID variable which stores, for each grid point of the ordered mesh, the corresponding node index of the unstructured dataset included in the input CFD solution. This allows mapping information to be used for different CFD solutions obtained on the same grid (e.g. for flow-field obtained at different time steps, in moving grid problems), with no need to re-compute nodes association.

You can use the following *BibTeX* citation as a reference (p. 318) for the **UNS2STR** tool:

@PhdThesis{MinervinoPhD,
  author   = {Minervino, Mauro},
  school   = {University of Naples Federico II},
  title    = {{A unified thermodynamic/vortical theory for the aerodynamic force analysis}},
  year     = {2025},
  month    = feb,
  type     = {phdthesis},
  doi      = {http://dx.doi.org/10.13140/RG.2.2.20705.83047},
  url      = {https://hdl.handle.net/11588/989058},
}

<ins>Contacts</ins>:
***Dr. Mauro Minervino***
*Established Researcher at the Italian Aerospace Research Centre (C.I.R.A. SCpA)*
m.minervino@cira.it
