

ratings <- qplot(rating, data=movies, geom="histogram")
b <- qplot(length, data=movies, geom="histogram")

png("a.png", width = output_width, height = output_height)
ggsave(file=filename)
dev.off()
