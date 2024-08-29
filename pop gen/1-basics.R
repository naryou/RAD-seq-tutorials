
# to get the working directory
getwd()

# to set the working directory
setwd('/Users/narcis/Narcis/Project/Courses_Conferences/Bioinformatics_NY/PartII/R_course/')

# to get your dataset in R
my_data <- read.table('./basics_ex.txt', header = TRUE, sep = '\t')

head(my_data) ########## To view the first rows
tail(my_data)

attach(my_data) ######To have direct access to all of the variables in the dataset

my_data[1,]
my_data[,1]
length(my_data)
dim(my_data)
nrow(my_data)
ncol(my_data)

# hist(my_data)  # on numeric columns
range(my_data [,2])  # on numeric columns
range(Length1)  ###Shows the range of data
sort(Length1)
length(Length1) #######Number of observations in dataset
summary(Length1) ####### Gives summary statistics
mean(Length1)
var(Length1)
sd(Length1)
sqrt(var(Length1)) ###The same as sd
hist(Length1, breaks=20, xlab = 'Length',col = "light blue", main='Histogram ... leaf') ########Histogram. R automatically decides about the breaks, but you can also specify it :)
###### main defines the title
??hist

round(sd(Length1), 2) ###To round the results, 1,2,... shows the decimal

boxplot(Length1, range = 0, main='Length', ylab='Length', border = rainbow(4)) #####without Range=0, the minimum and maximum interfaces (outliers?) are not shown
##### border defines the color of the boxplot

table(Length1) ######For categorical data, shows the number of each level
counts= table(Size)
sizes=c("L", "M", "S", "XL", "XS")
barplot(counts, col = rainbow(5), names.arg = sizes) ##### rainbow makes color of (here) 5
pie(counts, col=rainbow(5), label = sizes)

#
# Practice: basic math is easy in R
3+5                
12/7

weight_kg <- 55      # assigns values on the right to the object on the left

#Rules for objects
#Cannot start with a number; case sensitive; cannot be fundamental functions

2.2*weight_kg

weight_lb <- 2.2*weight_kg    
weight_kg <- 100                  #what is weight_lb now?

sample(1:41,1) #pick a number
#
