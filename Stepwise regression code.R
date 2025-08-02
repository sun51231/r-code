


library(clipr)
data <- read_clip_tbl()
data1<-data
colnames(data1) <- c("x1", "x2", "x3", "x4", "x5", "x6","x7", "y")
cor_data1<-cor(data1)
linear_model <- lm(y ~  x1  +x2  + x3 + x4+ x5 + x6 +x7, data = data1)
summary(linear_model)
step_model <- step(linear_model, direction = "both")  
summary(step_model)
