
dat <- dat_list[[90]]
incefreq2 <- buffers(dat, xy, nSite=1, r=150,
                     crs = 'epsg:4326', output = 'incidence_freq')

locs_inbuffer <- buffers(dat, xy, nSite = nSite, r = r, crs = crs, output = 'locs')
# Compute richness for each buffer using iNEXT::estimateD
# lapply is used to iterate over each buffer in incefreq2
richness_list <- lapply(incefreq2, function(buffer) {
  iNEXT::estimateD(x = buffer, q = c(0), datatype = "incidence_freq", 
                   base = "coverage", level = 0.7, nboot = 50)
})

# Combine the richness results from all iterations into a single data frame
# Iteration number is added as a column for tracking purposes
combined_richness <- do.call(rbind, lapply(seq_along(richness_list), function(i) {
  richness_list[[i]] %>%
    mutate(Iteration = i)  # Add iteration ID
}))

# Calculate summary statistics (mean and 95% confidence intervals) for each Assemblage
richness_summary <- combined_richness %>%
  group_by(Assemblage) %>%
  summarise(
    # Mean and confidence interval for `t`
    mean_t = mean(t),
    ci_t_lower = mean(t) - 1.96 * sd(t) / sqrt(n()),
    ci_t_upper = mean(t) + 1.96 * sd(t) / sqrt(n()),
    
    # Mean and confidence interval for `SC`
    mean_SC = mean(SC),
    ci_SC_lower = mean(SC) - 1.96 * sd(SC) / sqrt(n()),
    ci_SC_upper = mean(SC) + 1.96 * sd(SC) / sqrt(n()),
    
    # Mean and confidence interval for `qD`
    mean_qD = mean(qD),
    ci_qD_lower = mean(qD) - 1.96 * sd(qD) / sqrt(n()),
    ci_qD_upper = mean(qD) + 1.96 * sd(qD) / sqrt(n())
  )

# Apply Fun.ince to all buffers in all iterations of incefreq2
rich.info.all <- lapply(incefreq2, function(iteration) {
  # Apply Fun.ince to each buffer in the current iteration
  lapply(names(iteration), function(assembly_name) {
    result <- Fun.ince(iteration[[assembly_name]]) # Apply Fun.ince
    # Add Assemblage name as the first column
    c(Assemblage = assembly_name, result)
  })
})

# Combine results across all iterations into a single data frame
rich.infodf.all <- do.call(rbind, lapply(seq_along(rich.info.all), function(i) {
  iteration_results <- do.call(rbind, rich.info.all[[i]])
  # Add Iteration ID
  iteration_results <- cbind(Iteration = i, iteration_results)
  return(iteration_results)
}))

rich.infodf.all <- as.data.frame(rich.infodf.all)
# Set column names for the combined data frame
colnames(rich.infodf.all) <- c("Iteration", "Assemblage", "nT", "Number of occurrence",
                               "Sobs", "Chat", "Good.u", "multition.rate",
                               paste0("Q", 1:10))
# Columns to convert to numeric
numeric_cols <- c("nT", "Number of occurrence", "Sobs", "Chat", "Good.u", 
                  "multition.rate", paste0("Q", 1:10))

# Convert specified columns to numeric
rich.infodf.all[numeric_cols] <- lapply(rich.infodf.all[numeric_cols], as.numeric)
# Summarize results by buffer (Assemblage) across iterations
rich_info_summary <- rich.infodf.all %>%
  group_by(Assemblage) %>% 
  summarise(
    # Mean and confidence intervals for nT
    mean_nT = mean(nT),
    ci_nT_lower = mean(nT) - 1.96 * sd(nT) / sqrt(n()),
    ci_nT_upper = mean(nT) + 1.96 * sd(nT) / sqrt(n()),
    
    # Mean and confidence intervals for "Number of occurrence"
    mean_occurrence = mean(`Number of occurrence`),
    ci_occurrence_lower = mean(`Number of occurrence`) - 1.96 * sd(`Number of occurrence`) / sqrt(n()),
    ci_occurrence_upper = mean(`Number of occurrence`) + 1.96 * sd(`Number of occurrence`) / sqrt(n()),
    
    # Mean and confidence intervals for Sobs
    mean_Sobs = mean(Sobs),
    ci_Sobs_lower = mean(Sobs) - 1.96 * sd(Sobs) / sqrt(n()),
    ci_Sobs_upper = mean(Sobs) + 1.96 * sd(Sobs) / sqrt(n()),
    
    # Mean and confidence intervals for Chat
    mean_Chat = mean(Chat),
    ci_Chat_lower = mean(Chat) - 1.96 * sd(Chat) / sqrt(n()),
    ci_Chat_upper = mean(Chat) + 1.96 * sd(Chat) / sqrt(n()),
    
    # Mean and confidence intervals for Good.u
    mean_Good_u = mean(Good.u),
    ci_Good_u_lower = mean(Good.u) - 1.96 * sd(Good.u) / sqrt(n()),
    ci_Good_u_upper = mean(Good.u) + 1.96 * sd(Good.u) / sqrt(n()),
    
    # Mean and confidence intervals for multition.rate
    mean_multition_rate = mean(multition.rate),
    ci_multition_rate_lower = mean(multition.rate) - 1.96 * sd(multition.rate) / sqrt(n()),
    ci_multition_rate_upper = mean(multition.rate) + 1.96 * sd(multition.rate) / sqrt(n()),
    
    # Mean and confidence intervals for Q1-Q10
    across(starts_with("Q"), list(
      mean = ~ mean(.),
      ci_lower = ~ mean(.) - 1.96 * sd(.) / sqrt(n()),
      ci_upper = ~ mean(.) + 1.96 * sd(.) / sqrt(n())
    ), .names = "{.col}_{.fn}")
  )



# Manually create facet labels for Northern and Southern slopes
facet_labels <- data.frame(
  slope_type = c("Southern slope", "Northern slope"),
  x_label = 20,  # Rightmost position
  y_label = y_max_val * 0.85  # Position below max Y
)

# Create the slope variation curve plot
slope_vTime_plot <- ggplot() +
  # climate state color bars (Legend bar)
  geom_rect(data = climate_legend, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y*0.8, ymax=y*0.85),
            fill = climate_legend$climate_color) +  # Control bar size
  # Add "Cooler Climate" & "Warmer Climate" labels
  geom_text(aes(x = 355, y = y_max_val*0.8, label = "Cooler climate"), 
            size = 3, hjust = 1) +
  geom_text(aes(x = 270, y = y_max_val*0.8, label = "Warmer climate"), 
            size = 3, hjust = 0) +
  geom_rect(data = climate_states,
            aes(xmin = top, xmax = bottom, 
                ymin = y_max_val*0.9, ymax = y_max_val),
            fill = climate_states$climate_color,
            linewidth = 0.3) +
  # Add geological system color bands
  geom_rect(data = system_labels,
            aes(xmin = min_ma, xmax = max_ma, 
                ymin = y_min_val*0.9, ymax = y_min_val),
            fill = system_labels$fill_color,
            color = "black",
            linewidth = 0.3) +  
  # Add geological stage color bands
  geom_rect(data = time_bins,
            aes(xmin = min_ma, xmax = max_ma, 
                ymin = y_min_val*0.9, ymax = y_min_val*0.85),
            fill = time_bins$stageCol,
            color = "black",
            linewidth = 0.2) +
  # Remove fill color legend (to avoid redundancy)
  guides(fill = "none") +  
  # Add slope curves (different types)
  geom_line(data = slope_data, aes(x = bin_midpoint, y = slope_value, color = slope_type), linewidth = 1) +
  geom_point(data = slope_data, aes(x = bin_midpoint, y = slope_value, color = slope_type), size = 2) +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val), 
            color = "black", fill = NA, linewidth = 1) +
  geom_hline(yintercept = 0, color="black", linewidth=0.8) + 
  # Customize colors
  scale_color_manual(
    values = c("southern_slope" = "#E69F00", 
               "northern_slope" = "#0072B2"),
    labels = c("Southern slope", "Northern slope")
  ) +
  # Adjust x and y axis ranges
  scale_x_reverse(
    name = "Geological time (Ma)",
    limits = c(x_max_val, 0),  # Explicitly set limits to ensure alignment
    breaks = seq(500, 0, -50),
    expand = c(0, 0)  # Avoid extra whitespace
  ) +
  scale_y_continuous(limits = c(y_min_val, y_max_val),
                     expand = c(0, 0)                 # Disable axis expansion
  ) +
  # Add geological system names above the color bands
  geom_text(
    data = system_labels, 
    aes(x = mid_ma, y = y_text_pos, label = system),  # Use system names and computed midpoints
    angle = 0, 
    vjust = 0.5,      # Vertically centered
    hjust = 0.5,      
    size = 3
  ) +
  # Legends and labels
  labs(
    title = "Slopes of LDGs over geological time",
    x = "Geological time (Ma)",
    y = "Slope value",
    color = "Slope Type"
  ) +
  # Improve aesthetics using a clean theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "top",
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  ) 

# Save the plot as an image file
# ggsave("./figures/LDG slopes through geological time.jpg", slope_vTime_plot, width = 8, height = 5, dpi = 300)
# Display the plot
print(slope_vTime_plot)


### 3.4 vilion figure--------------------------------------------------------
# Convert data to long format and remove missing values (NA)
slope_melt <- slope_data %>%
  select(climate_state, slope_value, slope_type, everything()) %>%  # Ensure correct order
  rename(Slope_Value = slope_value, Slope_Type = slope_type) %>%  # Rename for consistency
  na.omit()  # Remove NA values


# Reorder the factor levels of climate_state
slope_melt$climate_state <- factor(slope_melt$climate_state, levels = climate_levels)

# Compute the total sample count for each climate state (combining all slope types)
sample_counts <- slope_cli_df_filter %>%
  select(climate_state, bin_midpoint) %>%
  group_by(climate_state) %>%  # Group only by climate state
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", count))  # Format label for display

# Adjust Y-axis position for text labels (using a fixed offset)
label_y_pos <- max(slope_melt$Slope_Value, na.rm = TRUE) + 0.01

# Create the violin plot with adjustments
violin_plot <- ggplot(slope_melt, aes(x = climate_state, y = Slope_Value, fill = Slope_Type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot without trimming
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +  # Overlay boxplot with no outliers
  # Add total sample count labels (centered above the violins)
  geom_text(
    data = sample_counts,
    aes(x = climate_state, y = label_y_pos, label = label),  # Use a fixed Y position
    inherit.aes = FALSE,  # Do not inherit fill mapping from Slope_Type
    size = 4,
    color = "black",
    vjust = -0.5  # Slight vertical adjustment
  ) +
  # Apply minimal theme
  theme_minimal() +
  labs(
    title = "Distribution of LDG slopes by climate state",
    x = "Climate state",
    y = "Slope calue",
    fill = "Slope type"
  ) +
  # Customize text and border styles
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 14),                             
    axis.text = element_text(size = 12),                              
    legend.title = element_text(size = 14),                         
    legend.text = element_text(size = 12),                         
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Save the violin plot as an image
# ggsave("./figures/LDG slope in climate states violin plot.jpg", violin_plot, width = 8, height = 5, dpi = 300)

# Display the plot
print(violin_plot)