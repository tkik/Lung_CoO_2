##################functions####################
label_func <- function(x){
  breaks <- x
  breaks[breaks==200] <- ">=200"
  breaks
}

region_annotation_function <- function(comparisons=comparisons_selected, gains=gain, losses=loss, annots_gr = annots_gr_orig){

  #browser()
plots <- list()
  for (comp in comparisons){

    factor_levels <- gsub("mm10_[[:alpha:]]*_", "", unique(annots_gr$type))
    factor_levels <- factor_levels[order(factor_levels, decreasing = T)]
    cat('\n')

    cat("### Comparisons", paste0(labels[labels$comparisons==comp,1], " ",
                                  labels[labels$comparisons==comp,2], " vs. ",
                                  labels[labels$comparisons==comp,3]), "\n")



    p_annots_data <- list()
    for (dataset in list(gains[[comp]], losses[[comp]])){


      direction_1 <- paste0("Higher in \n", labels[labels$comparisons==comp,1], " ", gsub("_", " ", labels[labels$comparisons==comp,3]),  " group")
      direction_2 <- paste0("Higher in \n", labels[labels$comparisons==comp,1], " ",
                            gsub("_", " ", labels[labels$comparisons==comp,2]),   " group")
      direction <- ifelse(all(mcols(dataset)$diff.Methy<0), direction_1, direction_2)

      result <- annotate_regions(dataset, annotations=annots_gr, minoverlap = 10L, ignore.strand = TRUE, quiet = FALSE)
      g <- getGenomeAndMask("mm10")$genome
      g <- filterChromosomes(g, chr.type="canonical", organism="mm10")
      rnd_annots = annotate_regions(regions = randomizeRegions(rep(dataset, 10), genome = g),annotations = annots_gr,ignore.strand = TRUE)

      p_annots = summarize_annotations(annotated_regions = result, annotated_random = rnd_annots)
      p_annots_data[[direction]] <- p_annots

      p_annots_data[[direction]]$annot.type <- gsub("[[:alnum:]]*_", "", p_annots_data[[direction]]$annot.type)

      p_annots_data[[direction]]$annot.type <-  factor(p_annots_data[[direction]]$annot.type, levels = factor_levels)

      p_annots_data[[direction]]$direction <- direction

    }


    dat <- do.call(rbind.data.frame, p_annots_data)
    dat <- dat %>% group_by(direction, data_type) %>% mutate(ratio=n/sum(n),sum=sum(n))


    # browser()
    ###########ggplot bars ############
    bar <- ggplot(dat)+
      geom_col(aes(x=annot.type, y = ratio, fill=data_type), position = "dodge")+ facet_grid(rows = vars(direction))+
      scale_fill_nejm()+
      theme_light()+
      theme(legend.position = "bottom", legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.8))+
      xlab("Annotation type")+
      scale_y_continuous(labels=scales::percent)+
      ylab("Percentage")+
      ggtitle(paste0("Distribution of regions with methylation alterations \n", labels[labels$comparisons==comp,1], " ",
                                  labels[labels$comparisons==comp,2], " vs. ",
                                  labels[labels$comparisons==comp,3]))

    print(bar)

    gain_regions <-  dat$annot.type[dat$direction==direction_1][duplicated(dat$annot.type[dat$direction==direction_1])]
    loss_regions <-   dat$annot.type[dat$direction==direction_2][duplicated(dat$annot.type[dat$direction==direction_2])]
    dat <- dat[(dat$direction==direction_1 & dat$annot.type %in% gain_regions) |
                 (dat$direction==direction_2 & dat$annot.type %in% loss_regions),]


    final <- dat %>% group_by(direction, annot.type) %>% summarise(ratio = log2(ratio[data_type=="Data"]/ratio[data_type=="Random Regions"]),
                                                                   chi_p = fisher.test(cbind(n, sum))$p.value)

    final$logp <- -log10(final$chi_p)
    final$logp[final$logp > 200] <- 200
    final$significant <- ifelse(final$chi_p<0.05, "yes", "no")
    final$significant <- factor(final$significant, levels=c("no", "yes"))

    #############ggplot_enrichment##################
    enr <- ggplot(data = final, aes(y=annot.type, x=direction))+
      geom_point(aes(size=logp, fill=ratio, color=significant), pch=21)+
      scale_fill_gradient2( midpoint = 0, low="#0072B5FF", high="#BC3C29FF", name = "Fold change")+
      scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
      scale_size(name="Log p", labels = label_func) +
      theme_bw()+
      theme(axis.text.x=element_text(size=12, angle = 90), axis.text.y=element_text(size=12),
            strip.text = element_text(size=12), aspect.ratio = 4)+
      ylab("Annotation type")+
      scale_x_discrete(name=NULL)+
      theme(legend.text=element_text(size=10),  legend.title=element_text(size=12))


    print(enr)

    cat("\n")
    plots[[comp]] <- list(bar, enr)
  }
return(plots)
}
