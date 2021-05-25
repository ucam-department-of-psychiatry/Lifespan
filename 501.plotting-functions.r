##
## Libraries (needed by functions in this file)
##
##

library("gamlss")

library("scales")

library("ggplot2")

##
## Plotting Functions
##
##


facet_zoom2 <- function(x, y, xy, zoom.data, xlim = NULL, ylim = NULL,
                        split = FALSE, horizontal = TRUE, zoom.size = 2,
                        show.area = TRUE, shrink = TRUE) {
    x <- if (missing(x)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(x)
    y <- if (missing(y)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(y)
    zoom.data <- if (missing(zoom.data)) NULL else lazyeval::lazy(zoom.data)
    if (is.null(x) && is.null(y) && is.null(xlim) && is.null(ylim)) {
        stop("Either x- or y-zoom must be given", call. = FALSE)
    }
    if (!is.null(xlim)) x <- NULL
    if (!is.null(ylim)) y <- NULL
    ggproto(NULL, FacetZoom2,
            shrink = shrink,
            params = list(
                x = x, y = y, xlim = xlim, ylim = ylim, split = split, zoom.data = zoom.data,
                zoom.size = zoom.size, show.area = show.area,
                horizontal = horizontal
            )
            )
}

## define FacetZoom as a ggproto object that inherits from FacetZoom,
## with a modified draw_panels function. the compute_layout function references
## the version currently on GH, which is slightly different from the CRAN
## package version.
FacetZoom2 <- ggproto(
    "FacetZoom2",
    ggforce::FacetZoom,
    compute_layout = function(data, params) {
        layout <- rbind( # has both x & y dimension
            data.frame(name = 'orig', SCALE_X = 1L, SCALE_Y = 1L),
            data.frame(name = 'x', SCALE_X = 2L, SCALE_Y = 1L),
            data.frame(name = 'y', SCALE_X = 1L, SCALE_Y = 2L),
            data.frame(name = 'full', SCALE_X = 2L, SCALE_Y = 2L),
            data.frame(name = 'orig_true', SCALE_X = 1L, SCALE_Y = 1L),
            data.frame(name = 'zoom_true', SCALE_X = 1L, SCALE_Y = 1L)
        )
        if (is.null(params$y) && is.null(params$ylim)) { # no y dimension
            layout <- layout[c(1,2, 5:6),]
        } else if (is.null(params$x) && is.null(params$xlim)) { # no x dimension
            layout <- layout[c(1,3, 5:6),]
        }
        layout$PANEL <- seq_len(nrow(layout))
        layout
    },
    draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                           data, theme, params) {
        if (is.null(params$x) && is.null(params$xlim)) {
            params$horizontal <- TRUE
        } else if (is.null(params$y) && is.null(params$ylim)) {
            params$horizontal <- FALSE
        }
        if (is.null(theme[['zoom']])) {
            theme$zoom <- theme$strip.background
        }
        if (is.null(theme$zoom.x)) {
            theme$zoom.x <- theme$zoom
        }
        if (is.null(theme$zoom.y)) {
            theme$zoom.y <- theme$zoom
        }
        axes <- render_axes(ranges, ranges, coord, theme, FALSE)
        panelGrobs <- ggforce:::create_panels(panels, axes$x, axes$y)
        panelGrobs <- panelGrobs[seq_len(length(panelGrobs) - 2)]
        if ('full' %in% layout$name && !params$split) {
            panelGrobs <- panelGrobs[c(1, 4)]
        }
                                        # changed coordinates in indicator / lines to zoom from
                                        # the opposite horizontal direction
        if ('y' %in% layout$name) {
            if (!inherits(theme$zoom.y, 'element_blank')) {
                zoom_prop <- scales::rescale(
                                         y_scales[[2]]$dimension(ggforce:::expansion(y_scales[[2]])),
                                         from = y_scales[[1]]$dimension(ggforce:::expansion(y_scales[[1]])))
                indicator <- polygonGrob(
                    x = c(0, 0, 1, 1), # was x = c(1, 1, 0, 0),
                    y = c(zoom_prop, 1, 0),
                    gp = gpar(col = NA, fill = alpha(theme$zoom.y$fill, 0.5)))
                lines <- segmentsGrob(
                    x0 = c(1, 1), x1 = c(0, 0), # was x0 = c(0, 0), x1 = c(1, 1)
                    y0 = c(0, 1), y1 = zoom_prop,
                    gp = gpar(col = theme$zoom.y$colour,
                              lty = theme$zoom.y$linetype,
                              lwd = theme$zoom.y$size,
                              lineend = 'round'))
                indicator_h <- grobTree(indicator, lines)
            } else {
                indicator_h <- zeroGrob()
            }
        }
        if ('x' %in% layout$name) {
            if (!inherits(theme$zoom.x, 'element_blank')) {
                zoom_prop <- scales::rescale(x_scales[[2]]$dimension(ggforce:::expansion(x_scales[[2]])),
                                             from = x_scales[[1]]$dimension(ggforce:::expansion(x_scales[[1]])))
                indicator <- polygonGrob(c(zoom_prop, 1, 0), c(1, 1, 0, 0),
                                         gp = gpar(col = NA, fill = alpha(theme$zoom.x$fill, 0.5)))
                lines <- segmentsGrob(x0 = c(0, 1), y0 = c(0, 0), x1 = zoom_prop, y1 = c(1, 1),
                                      gp = gpar(col = theme$zoom.x$colour,
                                                lty = theme$zoom.x$linetype,
                                                lwd = theme$zoom.x$size,
                                                lineend = 'round'))
                indicator_v <- grobTree(indicator, lines)
            } else {
                indicator_v <- zeroGrob()
            }
        }
        if ('full' %in% layout$name && params$split) {
            space.x <- theme$panel.spacing.x
            if (is.null(space.x)) space.x <- theme$panel.spacing
            space.x <- unit(5 * as.numeric(convertUnit(space.x, 'cm')), 'cm')
            space.y <- theme$panel.spacing.y
            if (is.null(space.y)) space.y <- theme$panel.spacing
            space.y <- unit(5 * as.numeric(convertUnit(space.y, 'cm')), 'cm')
                                        # change horizontal order of panels from [zoom, original] to [original, zoom]
                                        # final <- gtable::gtable_add_cols(panelGrobs[[3]], space.x)
                                        # final <- cbind(final, panelGrobs[[1]], size = 'first')
                                        # final_tmp <- gtable::gtable_add_cols(panelGrobs[[4]], space.x)
                                        # final_tmp <- cbind(final_tmp, panelGrobs[[2]], size = 'first')
            final <- gtable::gtable_add_cols(panelGrobs[[1]], space.x)
            final <- cbind(final, panelGrobs[[3]], size = 'first')
            final_tmp <- gtable::gtable_add_cols(panelGrobs[[2]], space.x)
            final_tmp <- cbind(final_tmp, panelGrobs[[4]], size = 'first')
            final <- gtable::gtable_add_rows(final, space.y)
            final <- rbind(final, final_tmp, size = 'first')
            final <- gtable::gtable_add_grob(final, list(indicator_h, indicator_h),
                                             c(2, 6), 3, c(2, 6), 5,
                                             z = -Inf, name = "zoom-indicator")
            final <- gtable::gtable_add_grob(final, list(indicator_v, indicator_v),
                                             3, c(2, 6), 5,
                                             z = -Inf, name = "zoom-indicator")
            heights <- unit.c(
                unit(max_height(list(axes$x[[1]]$top, axes$x[[3]]$top)), 'cm'),
                unit(1, 'null'),
                unit(max_height(list(axes$x[[1]]$bottom, axes$x[[3]]$bottom)), 'cm'),
                space.y,
                unit(max_height(list(axes$x[[2]]$top, axes$x[[4]]$top)), 'cm'),
                unit(params$zoom.size, 'null'),
                unit(max_height(list(axes$x[[2]]$bottom, axes$x[[4]]$bottom)), 'cm')
            )
                                        # swop panel width specifications according to the new horizontal order
            widths <- unit.c(
                                        # unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
                                        # unit(params$zoom.size, 'null'),
                                        # unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm'),
                                        # space.x,
                                        # unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
                                        # unit(1, 'null'),
                                        # unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')
                unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
                unit(1, 'null'),
                unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm'),
                space.x,
                unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
                unit(params$zoom.size, 'null'),
                unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm')
            )
            final$heights <- heights
            final$widths <- widths
        } else {
            if (params$horizontal) {
                space <- theme$panel.spacing.x
                if (is.null(space)) space <- theme$panel.spacing
                space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
                heights <- unit.c(
                    unit(max_height(list(axes$x[[1]]$top, axes$x[[2]]$top)), 'cm'),
                    unit(1, 'null'),
                    unit(max_height(list(axes$x[[1]]$bottom, axes$x[[2]]$bottom)), 'cm')
                )
                                        # change horizontal order of panels from [zoom, original] to [original, zoom]
                                        # first <- gtable::gtable_add_cols(panelGrobs[[2]], space)
                                        # first <- cbind(final, panelGrobs[[1]], size = 'first')
                final <- gtable::gtable_add_cols(panelGrobs[[1]], space)
                final <- cbind(final, panelGrobs[[2]], size = "first")
                final$heights <- heights
                                        # swop panel width specifications according to the new horizontal order
                                        # unit(c(params$zoom.size, 1), 'null')
                final$widths[panel_cols(final)$l] <- unit(c(1, params$zoom.size), 'null')
                final <- gtable::gtable_add_grob(final, indicator_h, 2, 3, 2, 5,
                                                 z = -Inf, name = "zoom-indicator")
            } else {
                space <- theme$panel.spacing.y
                if (is.null(space)) space <- theme$panel.spacing
                space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
                widths <- unit.c(
                    unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
                    unit(1, 'null'),
                    unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')
                )
                final <- gtable::gtable_add_rows(panelGrobs[[1]], space)
                final <- rbind(final, panelGrobs[[2]], size = 'first')
                final$widths <- widths
                final$heights[panel_rows(final)$t] <- unit(c(1, params$zoom.size), 'null')
                final <- gtable::gtable_add_grob(final, indicator_v, 3, 2, 5,
                                                 z = -Inf, name = "zoom-indicator")
            }
        }
        final
    }
)
