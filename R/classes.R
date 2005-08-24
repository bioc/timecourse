#	CLASSES.R

require(limma,quietly = FALSE, warn.conflicts = FALSE) 
require(MASS,quietly = FALSE, warn.conflicts = FALSE)
require(Biobase,quietly = FALSE, warn.conflicts = FALSE)
require(methods)

setClass("MArrayTC",
   representation("list")
)

setClass("LargeDataObject")
setIs("MArrayTC","LargeDataObject")


