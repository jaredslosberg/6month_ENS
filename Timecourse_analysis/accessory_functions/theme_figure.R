theme_figure <- function(){ 
  
 font <- ""   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      
      # #grid elements
      # panel.grid.major = element_blank(),    #strip major gridlines
      # panel.grid.minor = element_blank(),    #strip minor gridlines
      # axis.ticks = element_blank(),          #strip axis ticks
      
       #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 10,                #set font size
        face = 'bold',            #bold typeface
        vjust = 2),               #raise slightly

      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 8),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 8),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10))
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}