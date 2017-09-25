dist.fun <- function(x1,x2,y1,y2) {
  dist <- sqrt((y2-y1)^2 + (x2-x1)^2)
  return(dist)
}

randomwalk <- function(x) {
random.walk.x <- cumsum(runif(x))
random.walk.y <- cumsum(runif(x))

fl <- dist.fun(random.walk.x[1], random.walk.y[1], random.walk.x[x], random.walk.y[x])
distances <- c()
for(i in 1:(x-1)) {
  distances <- c(distances, dist.fun(random.walk.x[i], random.walk.y[i], random.walk.x[i+1], random.walk.y[i+1]))
}
distratio <- fl/sum(distances)
return(distratio)
}

walks <- replicate(9999, randomwalk(100))
hist(walks)
mean(walks)
sd(walks)
mean(walks) + 2*sd(walks)
