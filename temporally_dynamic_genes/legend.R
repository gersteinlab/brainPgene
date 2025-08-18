Monotonically_falling <- function(x){
  -2*x^2 + 1.28 * pi*x + log1p(10 - exp(1))
}

Monotonically_rising <- function(x){
  29.42 * x^2 + pi*x + log1p(10 - exp(1))
}

falling_then_rising <- function(x){
  0.39 * x^2 - pi*x + log1p(10 - exp(1))
}

rising_then_falling <- function(x){
  - 1.39 * x^2 + exp(1)*pi*x + log1p(10 - exp(1))
}

p1 <- ggplot(data.frame(x = c(9, 12)), aes(x = x)) +
  stat_function(fun = Monotonically_falling, col = pal_simpsons("springfield")(4)[1]) + labs(x = '', y = '') +
  theme_minimal() + 
  theme(axis.text=element_blank(), panel.grid = element_blank(),
        axis.ticks=element_blank(), panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black'))

p2 <- ggplot(data.frame(x = c(3.5, 6)), aes(x = x)) +
  stat_function(fun = Monotonically_rising, col = pal_simpsons("springfield")(4)[2]) + labs(x = '', y = '') +
  theme_minimal() + 
  theme(axis.text=element_blank(), panel.grid = element_blank(),
        axis.ticks=element_blank(), panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black'))

p3 <- ggplot(data.frame(x = c(3, 5.5)), aes(x = x)) +
  stat_function(fun = falling_then_rising, col = pal_simpsons("springfield")(4)[3]) + labs(x = '', y = '') +
  theme_minimal() + 
  theme(axis.text=element_blank(), panel.grid = element_blank(),
        axis.ticks=element_blank(), panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black'))

p4 <- ggplot(data.frame(x = c(2, 4)), aes(x = x)) +
  stat_function(fun = rising_then_falling, col = pal_simpsons("springfield")(4)[4]) + labs(x = '', y = '') +
  theme_minimal() + 
  theme(axis.text=element_blank(), panel.grid = element_blank(),
        axis.ticks=element_blank(), panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black'))

library(patchwork)
legend <- p1 + p2 + p3 + p4
ggsave('temporal_dynamic_expression/legend.pdf', legend, width = 5, height = 5)