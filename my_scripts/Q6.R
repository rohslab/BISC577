# Initialization
library(DNAshapeR)

# Extract sample sequences
fn <- "/Users/apple/BISC577HW2/bsun/hw2/CTCF/unbound_500.fa"

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
plotShape(pred$Roll)
plotShape(pred$HelT)
heatShape(pred$ProT, 20)


fn <- "/Users/apple/BISC577HW2/bsun/hw2/CTCF/bound_500.fa"

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
plotShape(pred$Roll)
plotShape(pred$HelT)
heatShape(pred$ProT, 20)
