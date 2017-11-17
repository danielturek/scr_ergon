Case study data used in SEPARATING MORTALITY AND EMIGRATION: MODELLING SPACE USE, DISPERSAL AND SURVIVAL WITH ROBUST-DESIGN SPATIAL-CAPTURE-RECAPTURE DATA by Torbjørn Ergon and Beth Gardner, 2013.

Data provided by Torbjørn Ergon, email: <t.h.ergon@ibv.uio.no>.

See Appendix S2 of the paper for detailed information about the data.

The file ‘ErgonAndGardner2013.rdat’ is an ASCII text representation of an R list object created by the ‘dput’ function in R. It can be read into R by the ‘dget’ function.

The data have the same structure as the data that can be simulated and fit to the RD-SCR model provided in Appendix S3 of the paper.

List elements are:

- N[1]: Number of individuals that are only seen in the last primary session or the session they are censored (these individuals are sorted first in the data)
- N[2]: Total number of individuals
- n.prim: Number of primary sessions
- dt[k]: Length of interval k in units of months (30 days)
- tod[k,j]: Time-of-day of trapping session k,j. 1 = evening, 2 = morning
- first[i]: First primary session of individual i
- K[i]: Last primary session for individual i (some individuals are censored)
- J[i,k]: Last secondary session for individual i in primary session k (some individuals are censored)
- gr[i]: Sex of individual i. 1 = female, 2 = male
- R: Number of traps
- X[r,]: Location of trap r = 1..R (x and y coordinate in meters)
- H[i,j,k]: Trap index for capture of individual i in secondary session j within primary session k. 1 = not captured, other values are r+1

‘dimnames’ attributes for the first dimension of ‘first’, ‘K’, ‘J’, ‘gr’ and ‘H’ are individual names in the full database.
