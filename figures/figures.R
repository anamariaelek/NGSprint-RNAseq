library(ggplot2)


# mapping
mapqc <- read.table("mapping-stats.tsv",fill=TRUE, header=TRUE)

ggplot(mapqc, aes(x=Aligner, y=Reads)) + 
  geom_boxplot() + 
  geom_point(size=2) + 
  geom_line(aes(group=Sample),alpha=0.2) + 
  scale_y_continuous(limits=c(60,100)) + 
  theme_light() + 
  labs(y="% mapped reads")
ggsave("Mapping.png",width=2.5,height=6)


ggplot(mapqc, aes(x=Aligner, y=Uniq_reads)) + 
  geom_boxplot() + 
  geom_point(size=2) + 
  geom_line(aes(group=Sample),alpha=0.2) + 
  scale_y_continuous(limits=c(60,100)) + 
  theme_light() + 
  labs(y="% uniquely mapped reads")
ggsave("Mapping-uniq.png",width=2.5,height=6)

# power
power <- read.table("power.tsv",fill=TRUE, header=TRUE)
ggplot(power, aes(Reads,DEGs)) + geom_point(size=2) + geom_smooth(method="lm") + scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6, accuracy=1))
ggsave("Power.png",width=2.5,height=2.5)
