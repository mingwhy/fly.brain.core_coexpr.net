#https://www.dataanalytics.org.uk/make-transparent-colors-in-r/#:~:text=Make%20transparent%20colors%20in%20R&text=The%20rgb()%20command%20is,255%20being%20%E2%80%9Csolid%E2%80%9D).

## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  #rgb.val <- col2rgb(color)
  
  rgb.val=col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  #invisible(t.col)
  t.col
}
## END
