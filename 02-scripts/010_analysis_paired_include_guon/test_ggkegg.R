library(ggkegg)
library(tidygraph)
library(dplyr)
graph <- ggkegg::pathway("hsa04110", use_cache=TRUE)
graph

graph |> 
  mutate(degree=centrality_degree(mode="all"),
         betweenness=centrality_betweenness()) |> 
  activate(nodes) |>
  filter(type=="gene") |>
  arrange(desc(degree)) |>
  as_tibble() |>
  relocate(degree, betweenness)


graph <- graph |> mutate(showname=strsplit(graphics_name, ",") |>
                           vapply("[", 1, FUN.VALUE="a"))

ggraph(graph, layout="manual", x=x, y=y)+
  geom_edge_parallel(aes(linetype=subtype_name),
                     arrow=arrow(length=unit(1,"mm"), type="closed"),
                     end_cap=circle(1,"cm"),
                     start_cap=circle(1,"cm"))+
  geom_node_rect(aes(fill=I(bgcolor),
                     filter=type == "gene"),
                 color="black")+
  geom_node_text(aes(label=showname,
                     filter=type == "gene"),
                 size=2)+
  theme_void()

graph |> mutate(x=NULL, y=NULL) |>
  ggraph(layout="nicely")+
  geom_edge_parallel(aes(color=subtype_name),
                     arrow=arrow(length=unit(1,"mm"), type="closed"),
                     end_cap=circle(0.1,"cm"),
                     start_cap=circle(0.1,"cm"))+
  geom_node_point(aes(filter=type == "gene"),
                  color="black")+
  geom_node_point(aes(filter=type == "group"),
                  color="tomato")+
  geom_node_text(aes(label=showname,
                     filter=type == "gene"),
                 size=3, repel=TRUE, bg.colour="white")+
  scale_edge_color_viridis(discrete=TRUE)+
  theme_void()

graph |>
  activate(nodes) |>
  mutate(hsa=convert_id("hsa")) |>
  filter(type == "gene") |>
  as_tibble() |>
  relocate(hsa)

graph |>
  activate(nodes) |>
  mutate(highlight=highlight_set_nodes("hsa:7157")) |>
  ggraph(layout="manual", x=x, y=y)+
  geom_node_rect(aes(fill=I(bgcolor),
                     filter=type == "gene"), color="black")+
  geom_node_rect(aes(fill="tomato", filter=highlight), color="black")+
  geom_node_text(aes(label=showname,
                     filter=type == "gene"), size=2)+
  geom_edge_parallel(aes(linetype=subtype_name),
                     arrow=arrow(length=unit(1,"mm"),
                                 type="closed"),
                     end_cap=circle(1,"cm"),
                     start_cap=circle(1,"cm"))+
  theme_void()

graph |>
  mutate(degree=centrality_degree(mode="all")) |>
  ggraph(graph, layout="manual", x=x, y=y)+
  geom_node_rect(aes(fill=degree,
                     filter=type == "gene"))+
  overlay_raw_map()+
  scale_fill_viridis_c()+
  theme_void()
