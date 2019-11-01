Class07
================
Zach Goldberg
10/23/2019

## Revisit functions fron last time

``` r
source("http://tinyurl.com/rescale-R")
```

Let’s try out the rescale
    function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale(c(3, 10, NA, 7))
```

    ## [1] 0.0000000 1.0000000        NA 0.5714286

Rescale2 has and stop()

``` r
rescale2(c(3, 10, NA, 7, "Zach"))
```

Write a function called both\_na that counts the number of positiions in
two vectors that are both NA

``` r
both_na <- function(x, y){
  sum(is.na(x) & is.na(y))
}

x <- c(NA, 2, 3, 4, NA, 6, 7)
y <- c(1, 2, 3, 4, NA, 6, 7)
both_na(x,y)
```

    ## [1] 1

Testing both\_na on more inputs

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

Deal with recycling

``` r
both_na2 <- function(x,y){
  if(length(x) != length(y)){
    stop("Input vectors are different lengths")
  }
  sum(is.na(x) & is.na(y))
}

both_na2(x,y2)
```

Write function to grade HW

``` r
# student 1
s1 <- c(100, 100, 100, 100, 100, 100, 100, 90) 

# student 2
s2 <- c(100, NA, 90, 90, 90, 90, 97, 80)

v <- sort(s2, decreasing = FALSE, na.last = FALSE)
which.min(s2)
```

    ## [1] 8

``` r
grade <- function(x){ #Takes a vector, removes the lowest score, and computes the average
 x[is.na(x)] <- 0 #Replace NAs with 0
 mean(x[-which.min(x)]) #Drop the lowest score and find the average
}

grade(s2)
```

    ## [1] 91

Write a function to grade a whole class

``` r
allgrade <- read.csv("http://tinyurl.com/gradeinput")
row.names(allgrade) <- allgrade[,1]
allgrade <- allgrade[,-1]

apply(allgrade, 1, grade)
```

    ##  student-1  student-2  student-3  student-4  student-5  student-6 
    ##      91.75      82.50      84.25      84.25      88.25      89.00 
    ##  student-7  student-8  student-9 student-10 student-11 student-12 
    ##      94.00      93.75      87.75      79.00      86.00      91.75 
    ## student-13 student-14 student-15 student-16 student-17 student-18 
    ##      92.25      87.75      78.75      89.50      88.00      94.50 
    ## student-19 student-20 
    ##      82.75      82.75
