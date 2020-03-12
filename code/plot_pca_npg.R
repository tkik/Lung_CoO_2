plot_pca_npg <- function (pca_res, m = NULL, col_anno = NULL, shape_anno = NULL,
          pc_x = "PC1", pc_y = "PC2", show_labels = FALSE)
{
  X <- Y <- color_me <- shape_me <- row_names <- NULL
  pc_vars <- pca_res$var_explained
  pca_res <- as.data.frame(pca_res$PC_matrix)
  pca_res$row_names <- rownames(pca_res)
  x_lab <- paste0(pc_x, " [", pc_vars[pc_x], " %]")
  y_lab <- paste0(pc_x, " [", pc_vars[pc_y], " %]")
  if (!is.null(col_anno) || !is.null(shape_anno)) {
    if (!is(object = m, class2 = "methrix")) {
      stop("Please provde methrix object while using col_anno or shape_anno")
    }
    pd <- as.data.frame(colData(m))
    pd <- pd[rownames(pca_res), , drop = FALSE]
    pca_res <- cbind(pca_res, pd)
  }
  if (!is.null(col_anno)) {
    col_anno_idx <- which(colnames(pca_res) == col_anno)
    if (length(col_anno_idx) == 0) {
      stop(paste0(col_anno, " not found in provided methrix object"))
    }
    else {
      colnames(pca_res)[col_anno_idx] <- "color_me"
    }
  }
  if (!is.null(shape_anno)) {
    shape_anno_idx <- which(colnames(pca_res) == shape_anno)
    if (length(shape_anno_idx) == 0) {
      stop(paste0(shape_anno, " not found in provided methrix object"))
    }
    else {
      colnames(pca_res)[shape_anno_idx] <- "shape_me"
    }
  }
  pc_x_idx <- which(colnames(pca_res) == pc_x)
  pc_y_idx <- which(colnames(pca_res) == pc_y)
  colnames(pca_res)[c(pc_x_idx, pc_y_idx)] <- c("X",
                                                "Y")
  if (all(c("color_me", "shape_me") %in% colnames(pca_res))) {
    pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, color = color_me,
                                         shape = shape_me, label = row_names)) + geom_point(size = 3) +
      xlab(pc_x) + ylab(pc_y) + labs(color = col_anno,
                                     shape = shape_anno) + scale_color_nejm()
  }
  else if ("color_me" %in% colnames(pca_res)) {
    pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, color = color_me,
                                         label = row_names)) + geom_point(size = 3) + xlab(pc_x) +
      ylab(pc_y) + labs(color = col_anno) + scale_color_nejm()
  }
  else if ("shape_me" %in% colnames(pca_res)) {
    pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, shape = shape_me,
                                         label = row_names)) + geom_point(size = 3) + xlab(pc_x) +
      ylab(pc_y) + labs(shape = shape_anno)
  }
  else {
    pca_gg <- ggplot(data = as.data.frame(pca_res), aes(x = X,
                                                        y = Y, label = row_names)) + geom_point(size = 3,
                                                                                                fill = "black", color = "gray70") + xlab(pc_x) +
      ylab(pc_y)
  }
  pca_gg <- pca_gg + xlab(label = x_lab) + ylab(label = y_lab) +
    theme_classic(base_size = 12) + theme(axis.text.x = element_text(colour = "black",
                                                                     size = 12), axis.text.y = element_text(colour = "black",
                                                                                                            size = 12))
  if (show_labels) {
    pca_gg <- pca_gg + geom_label(size = 4)
  }
  pca_gg
}
