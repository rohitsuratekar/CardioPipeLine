library("ggplot2")
d <-  data.frame(price= runif(10, min=10, max=15))
g <- ggplot(d, aes(x=factor(1:length(d$price)) , y=price)) +
  geom_bar(stat="identity", fill="steelblue") +
  ggtitle("Example") +
  xlab("Index") +
  expand_limits(y=c(0,20)) +
  scale_y_continuous(expand = c(0,0))+
  theme_bw() +
  theme(panel.border = element_rect(fill=NA), 
        panel.background = element_rect(fill = alpha("gray", 0.6)))
plot(g)