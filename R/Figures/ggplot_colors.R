cvi_colours = list(
  cvi_std = c("#03071e",
    "#370617", "#6a040f",
    "#9d0208", "#d00000",
    "#dc2f02", "#e85d04",
    "#f48c06", "#faa307",
    "#ffba08"),
  blue_red = c("#b7094c", "#a01a58",
             "#892b64", "#723c70",
             "#5c4d7d", "#455e89",
             "#2e6f95", "#1780a1","#0091ad"),
  my_favourite_colours = c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51"),
  halloween = c("#9836e8", "#50db6e", "#f36f32", "#1b1f1f")
)
cvi_palettes = function(name, n, all_palettes = cvi_colours, type = c("discrete", "continuous")) {
  palette = all_palettes[[name]]
  if (missing(n)) {
    n = length(palette)
  }
  type = match.arg(type)
  out = switch(type,
               continuous = grDevices::colorRampPalette(palette)(n),
               discrete = palette[1:n]
  )
  structure(out, name = name, class = "palette")
}
scale_colour_cvi_d = function(name) {
  ggplot2::scale_colour_manual(values = cvi_palettes(name,
                                                     type = "discrete"))
}

scale_fill_cvi_d = function(name) {
  ggplot2::scale_fill_manual(values = cvi_palettes(name,
                                                   type = "discrete"))
}

scale_colour_cvi_c = function(name) {
  ggplot2::scale_colour_gradientn(colours = cvi_palettes(name = name,
                                                         type = "continuous"))
}

scale_fill_cvi_c = function(name) {
  ggplot2::scale_fill_gradientn(colours = cvi_palettes(name = name,
                                                       type = "continuous"))
}
