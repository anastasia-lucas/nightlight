#' gxeman
#'
#' Create heatmap plots for SNPxE or GxE interaction
#' Dependencies: ggplot2
#' @param d data frame, must contain SNP, ENV, CHR, POS, pvalue columns, optional Group column (for grouping ENV)
#' @param snpme optional main effect data frame containing SNP, CHR, POS, pvalue columns, optional Shape column
#' @param envme optional main effect data frame containing ENV and pvalue columns, optional Group column (for grouping ENV)
#' @param highlight_p threshold to diverge gradient color scale, default 0.05, set to "off" for no threshold
#' @param legend title for color legend, default "p-value"
#' @param title optional string for plot title
#' @param high color for high values
#' @param low color for low values
#' @param highlight_high if highlight_p given, color for max of highlight range
#' @param highlight_low if highlight_p given, color for min of highlight range
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image
#' @export
#' @examples
#' gxeman(d, me, symmetric, highlight_p, legend, title, high, low, highlight_high, highlight_low, file, hgt, wi, res)

gxeman <- function(d, snpme, envme, highlight_p=0.05, legend="p-value", title=NULL, high="#02021e", low="#1A0F99", highlight_high="yellow", highlight_low="#fffcd3", file="gxgman", hgt=7, wi=7.5, res=300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 and ggiraph to create interactive visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly=TRUE)
  }

  pvals <- d$pvalue
  if(!missing(snpme)){
    pvals <- c(pvals, envme$pvalue)
  }
  if(!missing(envme)){
    pvals <- c(pvals, snpme$pvalue)
  }

  d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))

  #Initialize vector for gradient
  if(highlight_p!="off"){
    color_vector <- c(highlight_low, highlight_high, low, high)
    #Values should be min, max < highlight value, highlight value, and max
    value_vector <- c(floor(min(pvals)), max(pvals[pvals < highlight_p]), highlight_p,ceiling(max(pvals)))
    breaks_vector <- c(floor(min(pvals)), max(pvals[pvals < highlight_p]), highlight_p, ceiling(max(pvals)))
    label_vector <- c(floor(min(pvals)), paste("<",highlight_p,sep=""), highlight_p, ceiling(max(pvals)))
  } else {
    color_vector <- c(low, high)
    #Values should be min, max < highlight value, highlight value, and max
    value_vector <- c(floor(min(pvals)), ceiling(max(pvals)))
    breaks_vector <- c(floor(min(pvals)), ceiling(max(pvals)))
    label_vector <- c(floor(min(pvals)), ceiling(max(pvals)))
  }

  #Create position index for chr:bp
  dsub <- d[, colnames(d) %in% c("SNP", "CHR", "POS")]
  dsuborder <- unique(dsub[order(dsub$CHR, dsub$POS), ])
  dsuborder$snp_pos_index <- seq.int(nrow(dsuborder))
  dpos <- merge(dsuborder[, colnames(dsuborder) %in% c("SNP", "snp_pos_index")], d, by="SNP", all.y=TRUE)

  #Create position index for non-genetic data
  esub <- d[, colnames(d) %in% c("ENV", "Group")]
  if("Group" %in% colnames(esub)){
    esub$ENVcol <- esub$Group
  } else {
    esub$ENVcol <- esub$ENV
  }
  esuborder <- unique(esub[order(esub$ENVcol), ])
  esuborder$env_pos_index <- seq.int(nrow(esuborder))
  epos <- merge(esuborder[, colnames(esuborder) %in% c("ENV", "ENVcol", "env_pos_index")], dpos, by="ENV", all.y=TRUE)

  #Set up axis for chr:bp
  xmaxRows <- by(dsuborder, dsuborder$CHR, function(x) x[which.max(x$snp_pos_index),])
  xminRows <- by(dsuborder, dsuborder$CHR, function(x) x[which.min(x$snp_pos_index),])
  xmilimits <- do.call(rbind, xminRows)
  xmalimits <- do.call(rbind, xmaxRows)
  xlims <- merge(xmilimits, xmalimits, by="CHR")
  names(xlims) <- c("CHR", "snpx", "posx", "posmin", "snpy", "posy", "posmax")
  xlims$av <- (xlims$posmin + xlims$posmax)/2
  xlims <- xlims[order(xlims$CHR),]

  #Set up axis for non-genetic data
  ymaxRows <- by(esuborder, esuborder$ENVcol, function(x) x[which.max(x$env_pos_index),])
  yminRows <- by(esuborder, esuborder$ENVcol, function(x) x[which.min(x$env_pos_index),])
  ymilimits <- do.call(rbind, yminRows)
  ymalimits <- do.call(rbind, ymaxRows)
  ylims <- merge(ymilimits[, -which(names(ymilimits)=="Group")], ymalimits[, -which(names(ymalimits)=="Group")], by="ENVcol")
  names(ylims) <- c("ENVcol", "envx", "posmin", "envy", "posmax")
  ylims$av <- (ylims$posmin + ylims$posmax)/2
  ylims <- ylims[order(ylims$ENVcol),]

  #Plot
  p <- ggplot(epos, aes(x=snp_pos_index, y=env_pos_index, fill=pvalue)) + geom_tile()
  p <- p + scale_x_continuous(breaks=xlims$posmax+0.5, labels=xlims$CHR, expand=c(0,0))
  p <- p + scale_y_continuous(breaks=ylims$posmax+0.5, labels=ylims$ENVcol, expand=c(0,0))
  #Optional add SNP main effect
  if(!missing(snpme) | !missing(envme)){
    if(!missing(snpme)){
      #Get position
      xme <- merge(dsuborder[, colnames(dsuborder) %in% c("SNP", "snp_pos_index")], snpme, by="SNP", all=TRUE)
      if("Shape" %in% names(snpme)){
        xme$Shape <- factor(xme$Shape)
        p <- p + geom_point(data=xme, aes(x=snp_pos_index, y= -1, color=pvalue, shape=Shape)) + expand_limits(y=-2)
        p <- p + labs(shape="")
      } else {
        p <- p + geom_point(data=xme, aes(x=snp_pos_index, y= -1, color=pvalue)) + expand_limits(y=-2)
      }
    }
    if(!missing(envme)){
      #Get position
      yme <- merge(esuborder[, colnames(esuborder) %in% c("ENV", "ENVcol", "env_pos_index")], snpme, by="ENV", all=TRUE)
      if("Shape" %in% names(snpme)){
        yme$Shape <- factor(yme$Shape)
        p <- p + geom_point(data=yme, aes(x=-1, y=env_pos_index, color=pvalue, shape=Shape)) + expand_limits(x=-2)
        p <- p + labs(shape="")
      } else {
        p <- p + geom_point(data=yme, aes(x=-1, y=env_pos_index, color=pvalue)) + expand_limits(x=-2)
      }
    }
    #Fix guides
    p <- p + scale_color_gradientn(colours=color_vector, values=value_vector, name=legend,limits=c(0,1),breaks=breaks_vector, labels=label_vector, FALSE)
    p <- p + guides(fill = guide_legend(override.aes = list(shape = NA)))
  }
  #Add colorbar or legend
  if(highlight_p!="off"){
    p <- p + scale_fill_gradientn(colours=color_vector, values=value_vector, name=legend,limits=c(0,1),breaks=breaks_vector, labels=label_vector, guide="legend")
  } else {
    p <- p + scale_fill_gradientn(colours=color_vector, values=value_vector, name=legend,limits=c(0,1),breaks=breaks_vector, labels=label_vector, guide="colorbar")
  }
  #Format
  p <- p + theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  p <- p + ggtitle(title)

  return(p)
}
