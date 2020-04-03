# Reproducible example taken from skbio's Unifract function
#   http://scikit-bio.org/docs/0.4.2/generated/generated/skbio.diversity.beta.unweighted_unifrac.html

library(phyloseq)
library(ape)

# Make OTU table
u_counts = c(1, 0, 0, 4, 1, 2, 3, 0)
v_counts = c(0, 1, 1, 6, 0, 1, 0, 0)
table=data.frame(u=u_counts, v=v_counts, row.names=paste("OTU", 1:8, sep=""))
tree = read.tree(text="(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,(OTU5:0.5,((OTU6:0.33,OTU7:0.62):0.5,OTU8:0.5):0.5):0.5):1.25):0.0)root;")
mydata = phyloseq(otu_table(table, taxa_are_rows=TRUE), phy_tree(tree))

# Unifrac; unweighted should be 0.37
UniFrac(mydata, weighted=FALSE)
