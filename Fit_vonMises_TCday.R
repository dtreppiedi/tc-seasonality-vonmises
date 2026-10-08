######################################################################
####   FIT OF A VON MISES MIXTURE TO A TCday SAMPLE               ####
######################################################################
#
#  Treppiedi, D., Villarini, G., Scoccimarro, E.,
#  Yang, W., Vecchi, G. A., Noto, L. V. - Contrasting Shifts in Tropical
#  Cyclone Seasonality Across Global Basins, Geophysical Research Letters.
#
#  It reads one sample of TCday (one basin, one intensity class, one time
#  window) from a csv file, fits a mixture of von Mises distributions and
#  returns the percentiles used in the paper (TCnth) and the figure.
#
#  The day of the year is a circular variable, so the fit is done on the
#  circle, with one or two components chosen by the BIC. The circle is cut
#  at the day of minimum density of the basin, so that the 5th percentile
#  is the onset of the season and the 95th is the decay.
#
######################################################################

rm(list = ls(all=TRUE))

library("movMF")
library("ggplot2")
library("stringr")


####      SETTINGS      ####

# input file, one row per TCday; only the column Date is used

NAME_CSV = "TCdays_NorthAtlantic_2010_2024.csv"

# where the circle is cut, i.e. the day with the lowest density of TCday.

DAY_MIN_DENS = "2001-03-01"

PERCENTILES = c(0.05, 0.25, 0.50, 0.75, 0.95)

NMIX_MAX = 2          # components tried in the mixture, the best by BIC

NSAMPLE = 100000      # draws used to invert the fitted distribution


####      EXTERNAL FUNCTIONS      ####

##------------------------------COMPUTE FRACTION OF YEAR------------##
# return fraction of year;
# argument TS, a dataframe with a "Date" column

year.fraction <- function(TS){
  
  ind_date <- as.Date(TS$Date, tz = "GMT")
  
  year_date = as.numeric(format(ind_date, "%Y"))
  
  is_leap = (year_date %% 4 == 0 & year_date %% 100 != 0) | (year_date %% 400 == 0)
  
  days_in_year = ifelse(is_leap, 366, 365)
  
  start_years <- as.Date(paste(year_date, "01", "01", sep = "-"), tz = "GMT")
  
  ind_jdate <- as.numeric(ind_date - start_years)
  
  ind_jdate = ind_jdate/days_in_year
  
  return(ind_jdate)
  
}


##-----------------FIT MIXTURE OF VON MISES DISTRIBUTIONS------------##
# function to fit a mixture of von Mises distributions
# returns a fitted movMF object
# argument fraclist = list of year fractions
# argument nmix = number of models to fit in the mixture

fit.mix.vMF.mod <- function(fraclist, nmix = 1){

  angles <- 2*pi*fraclist[!is.na(fraclist)]  # remove NAs to keep movMF happy
  hy = sin(angles)  ;  hx = cos(angles)

  set.seed(42)  # the answer to everything and to ensure reproducible results

  mises <- tryCatch(movMF(cbind(hx,hy), nmix, control = list(nruns = 20)), error=function(e){
    if (length(e)>0) print(e); return(NA)})

}


#------------Density of Fitted Mixture of VMFs--------------#
# argument is a movMF object

vMF.density.mod <- function(vMF){

  angles <- 2*pi*seq(0,1,length.out = 720)  # angles for every half degree
  xy <- cbind(cos(angles), sin(angles))

  return(data.frame(angles, density=dmovMF(xy, vMF$theta, vMF$alpha)))

}


####      1. READ THE TCday SAMPLE      ####

DATA = read.csv(NAME_CSV, stringsAsFactors = FALSE)

DATA$Date = as.Date(DATA$Date)

DATA$Fraction = year.fraction(DATA)

DATA$Angles = DATA$Fraction * 2 * pi

print(paste("TCday read:", nrow(DATA),
            "- from", min(DATA$Date), "to", max(DATA$Date)))


a_min = year.fraction(data.frame(Date = as.Date(DAY_MIN_DENS))) * 2 * pi
a_max = a_min + 2 * pi

DATA$Angles_shifted = DATA$Angles

DATA$Angles_shifted[which(DATA$Angles < a_min)] = DATA$Angles[which(DATA$Angles < a_min)] + 2 * pi


####      2. FIT THE MIXTURE AND CHOOSE THE NUMBER OF COMPONENTS      ####

VM_LIST = list()
BIC_ALL = rep(NA, NMIX_MAX)

for (kk in 1:NMIX_MAX) {

  VM_LIST[[kk]] = fit.mix.vMF.mod(DATA$Fraction, nmix = kk)

  if(inherits(VM_LIST[[kk]], "movMF")){

    BIC_ALL[kk] = stats::BIC(VM_LIST[[kk]])

    print(paste("components:", kk, "- BIC:", round(BIC_ALL[kk], 1)))

  }else{

    BIC_ALL[kk] = Inf

    print(paste("components:", kk, "- the fit did not converge"))

  }

}

if(all(is.infinite(BIC_ALL))) stop("no mixture could be fitted to this sample")

best_mix = which.min(BIC_ALL)

VM_BEST = VM_LIST[[best_mix]]

print(paste("best mixture:", best_mix, "component(s)"))


# density of the fitted mixture, to be drawn on the circle

toplotdens = vMF.density.mod(VM_BEST)


# mean direction and weight of each component, the arrows in the figure

toplotvect = data.frame(Pos = ifelse((atan2(VM_BEST$theta[,2], VM_BEST$theta[,1])/(2*pi)) < 0,
                                     yes = 1 + (atan2(VM_BEST$theta[,2], VM_BEST$theta[,1])/(2*pi)),
                                     no  = (atan2(VM_BEST$theta[,2], VM_BEST$theta[,1])/(2*pi))),
                        Len = VM_BEST$alpha)


####      3. PERCENTILES OF THE FITTED DISTRIBUTION      ####

set.seed(42)

# 1) Generate 100k random angles from the best fitted vMF model

temp = rmovMF(NSAMPLE, theta = VM_BEST$theta, alpha = VM_BEST$alpha)
temp = data.frame(x = temp[,1], y = temp[,2])
temp$ang = atan2(temp$y, temp$x)

# 2) Generate 2 identical sample and add 2*pi and 4*pi; this allow to linearize
#    a circular sample, it is as if I unwrap 3 times the circular sample,
#    keeping the continuity. It is necessary because rmovMF generates points
#    between -pi and pi

temp_forw1 = temp
temp_forw2 = temp

temp_forw1$ang = temp_forw1$ang + 2*pi
temp_forw2$ang = temp_forw2$ang + 4*pi

temp_forw_copang = rbind(temp, temp_forw1, temp_forw2)

# 3) I keep the generated angles between the cut of the circle (a_min) and
#    that angle plus 2*pi (a_max)

vMsCDF_forw_copang = temp_forw_copang[which(temp_forw_copang$ang >= a_min &
                                              temp_forw_copang$ang <  a_max),]

# 4) now 'quantile' gives the percentiles. It is the empirical quantile of a
#    100k sample generated from the distribution, so it is the theoretical one

TEORQUAN = quantile(vMsCDF_forw_copang$ang, PERCENTILES)

TEORQUAN.df = data.frame(Percentile = names(TEORQUAN),
                         Angle      = as.numeric(TEORQUAN),
                         Fraction   = as.numeric(TEORQUAN)/(2*pi))

# the same percentile written as a date and as the day of the season

TEORQUAN.df$Date = format(as.Date("1979-01-01") +
                            round((TEORQUAN.df$Fraction %% 1) * 365), "%d %b")

TEORQUAN.df$Day_of_season = round((TEORQUAN.df$Angle - a_min) / (2*pi) * 365, 1)

# height of each percentile on the density curve, it serves for the figure

TEORQUAN.df$density = NA

for (ii in 1:nrow(TEORQUAN.df)) {

  pos = which.min(abs(toplotdens$angles - (TEORQUAN.df$Fraction[ii] %% 1) * 2 * pi))

  TEORQUAN.df$density[ii] = toplotdens$density[pos]

}

print(TEORQUAN.df[, c("Percentile", "Date", "Day_of_season")])


####      4. FIGURE      ####

fract_label = data.frame(Fract = year.fraction(data.frame(Date = as.Date(str_c(month.abb, "01", "1979", sep = "-"),
                                                                        format = "%b-%d-%Y"))),
                         y = NA)

ymax = 2 * ceiling(max(toplotdens$density))

fract_label$y = -0.18 * ymax   # month labels, just inside the circle

PLOT_SEASON = ggplot()+
  geom_histogram(DATA, mapping = aes(x = Angles/(2*pi), y = after_stat(density)), boundary = 0,
                 fill = "gray30", alpha = 0.3, binwidth = 1/365)+
  geom_line(toplotdens, mapping = aes(x = angles/(2*pi), y = density),
            linewidth = 1, color = "red")+
  geom_segment(TEORQUAN.df, mapping = aes(x = Fraction %% 1, xend = Fraction %% 1,
                                          y = 0, yend = density, color = Percentile),
               linewidth = 1)+
  geom_point(TEORQUAN.df, mapping = aes(x = Fraction %% 1, y = density, fill = Percentile),
             size = 3, shape = 21)+
  geom_segment(toplotvect, mapping = aes(x = Pos, xend = Pos,
                                         y = -ymax, yend = -ymax + Len*ymax),
               linewidth = 1, color = "red",
               arrow = arrow(length = unit(0.2, "cm")))+
  scale_fill_manual(name = "Percentile",
                    breaks = c("5%", "25%", "50%", "75%", "95%"),
                    values = c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4"),
                    drop = FALSE)+
  scale_color_manual(name = "Percentile",
                     breaks = c("5%", "25%", "50%", "75%", "95%"),
                     values = c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4"),
                     drop = FALSE)+
  geom_text(fract_label, mapping = aes(x = Fract, y = y, label = month.abb), size = 4)+
  geom_segment(fract_label, mapping = aes(x = Fract, xend = Fract, y = 0, yend = -0.05*ymax))+
  geom_hline(yintercept = 0, color = "black")+
  xlim(0, 1)+
  ylim(-ymax, ymax)+
  coord_polar(start = 0)+
  theme_minimal()+
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

print(PLOT_SEASON)


####      5. SAVE      ####

write.csv(TEORQUAN.df[, c("Percentile", "Date", "Day_of_season", "Fraction")],
          "Percentiles.csv", row.names = FALSE)

ggsave(filename = "Seasonality.png", plot = PLOT_SEASON,
       width = 6.5, height = 5, dpi = 300)

print("done")
