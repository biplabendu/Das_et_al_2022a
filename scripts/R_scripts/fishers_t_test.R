a <- c(1:100)
b <- c(51:250)

# making the contingency table for a Fisgher's exact test
# find length of overlap
x1 <- intersect(a,b) %>% length()
# find numbers of genes in a that are not in b
x2 <- setdiff(a,b) %>% length()

# not in a but in b
x3 <- setdiff(b,a) %>% length()

# not in a and not in b
x4 <- 8000

c.tab <- data.frame(c1 = c(x1,x2),
                    c2 = c(x3,x4))
c.tab
fisher.test(c.tab)
