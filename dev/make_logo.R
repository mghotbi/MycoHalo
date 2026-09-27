W <- 1037; H <- 1200
png("logo.png", W, H, bg = "transparent", type = "cairo", antialias = "subpixel")
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot.new(); plot.window(c(-1, 1), c(-1.1575, 1.1575), asp = 1)
hex <- function(r) { a <- pi/2 + (0:5) * pi/3; list(x = r * cos(a), y = r * sin(a)) }
polygon(hex(1.1575), col = "#C9971C", border = NA)
polygon(hex(1.1575 - 0.07), col = "#2B1D16", border = NA)
cx <- 0; cy <- 0.16; R <- 0.64
for (k in 40:1) symbols(cx, cy, circles = 0.3 + k * 0.011, inches = FALSE, add = TRUE,
                        bg = adjustcolor("#6D4C41", alpha.f = 0.035), fg = NA)
circ <- function(x, y, r, n = 300) { t <- seq(0, 2*pi, length.out = n); list(x = x + r*cos(t), y = y + r*sin(t)) }
polygon(circ(cx, cy - 0.025, R + 0.04), col = adjustcolor("black", 0.35), border = NA)
polygon(circ(cx, cy, R + 0.04), col = "#D8DCE0", border = NA)
polygon(circ(cx, cy, R), col = "#3A3F45", border = NA)
for (k in 30:1) polygon(circ(cx, cy, 0.13 + k * 0.0085), col = adjustcolor("#E9B949", alpha.f = 0.045), border = NA)
rfun <- function(r, amp, seed) { set.seed(seed); k <- 3:8
  a <- runif(length(k), 0, amp); p <- runif(length(k), 0, 2*pi)
  function(t) r * (1 + Reduce(`+`, lapply(seq_along(k), function(j) a[j]*cos(k[j]*t + p[j])))) }
blob <- function(x, y, rf, n = 400) { t <- seq(0, 2*pi, length.out = n); rr <- rf(t); list(x = x + rr*cos(t), y = y + rr*sin(t)) }
rb <- rfun(0.12, 0.05, 3)
polygon(blob(cx, cy, rb), col = "#FFF3D6", border = "#E8D3A0", lwd = 3)
# colonies as smooth rasters: grey, melanized (sepia-dark) towards the bacterium
pos <- list(c(-1, 1), c(1, 1), c(1, -1), c(-1, -1)); d <- 0.30
base <- rbind(c(160,156,150), c(146,140,133), c(132,124,116), c(170,165,158)) / 255
dark <- c(78, 52, 46) / 255
face <- function(x, y) {
  points(c(x - 0.04, x + 0.04), c(y + 0.012, y + 0.012), pch = 19, cex = 3, col = "#1B1210")
  points(c(x - 0.032, x + 0.048), c(y + 0.024, y + 0.024), pch = 19, cex = 1, col = "white")
  t <- seq(pi*1.15, pi*1.85, length.out = 30)
  lines(x + 0.035*cos(t), y - 0.012 + 0.035*sin(t), lwd = 5, col = "#1B1210")
  points(c(x - 0.085, x + 0.085), c(y - 0.03, y - 0.03), pch = 19, cex = 2.6, col = adjustcolor("#F48FB1", 0.75)) }
set.seed(7)
for (i in 1:4) {
  x <- cx + pos[[i]][1] * d; y <- cy + pos[[i]][2] * d
  rf <- rfun(0.15, 0.07, 10 + i)
  n <- 260; gx <- seq(x - 0.2, x + 0.2, length.out = n); gy <- seq(y + 0.2, y - 0.2, length.out = n)
  X <- matrix(gx, n, n, byrow = TRUE); Y <- matrix(gy, n, n)
  th <- atan2(Y - y, X - x); rr <- sqrt((X - x)^2 + (Y - y)^2)
  inside <- rr <= rf(th)
  ang <- atan2(cy - y, cx - x)
  proj <- ((X - x) * cos(ang) + (Y - y) * sin(ang)) / 0.15
  w <- plogis((proj - 0.25) * 7) * 0.9
  w <- w + 0.08 * (1 - pmin(1, rr / rf(th)))           # darker centre (older tissue)
  tex <- matrix(rnorm(n * n), n, n); tex <- (tex + tex[c(2:n, 1), ] + tex[, c(2:n, 1)]) / 3
  w <- pmin(1, pmax(0, w + 0.05 * tex))
  edge <- pmin(1, (rf(th) - rr) / 0.02); a <- ifelse(inside, pmin(1, edge), 0)
  rgbm <- lapply(1:3, function(k) base[i, k] * (1 - w) + dark[k] * w)
  img <- matrix(rgb(rgbm[[1]], rgbm[[2]], rgbm[[3]], pmax(a, 0)), n, n)
  rasterImage(as.raster(img), x - 0.2, y - 0.2, x + 0.2, y + 0.2, interpolate = TRUE)
  polygon(blob(x, y, rf), col = NA, border = adjustcolor("#F5EFE6", 0.6), lwd = 3)
  face(x - 0.02 * cos(ang), y - 0.02 * sin(ang))
}
spark <- function(x, y, s) { segments(x - s, y, x + s, y, col = "#FFE9A8", lwd = 5); segments(x, y - s, x, y + s, col = "#FFE9A8", lwd = 5) }
spark(0.6, 0.66, 0.03); spark(-0.68, -0.33, 0.025); spark(0.66, -0.38, 0.02)
text(0, -0.72, "MycoHalo", col = "#F6E7C8", cex = 8.2, font = 2)
dev.off()
