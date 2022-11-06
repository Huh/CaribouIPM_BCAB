## ----render, eval=FALSE,include=FALSE----------------------------------------------------------------------------------------------------------------------
## rmarkdown::render(here::here('CaribouIPM_BCAB.Rmd'),
##                   output_file =  "README.md")
## 
## knitr::purl(input=here::here('CaribouIPM_BCAB.Rmd'),
##                   output =  here::here("scripts_r", '3.results_summary.r'))


## ----Load packages and data, results='hide', message=FALSE, warning=FALSE----------------------------------------------------------------------------------
library(ggmap)
library(raster)
library(RStoolbox)
library(ggsn)
library(MCMCvis)
library(tidybayes)
library(ggmcmc)
library(boot)
library(cowplot)
library(hrbrthemes)
library(lme4)
library(tidymodels)
library(broom.mixed)
library(RColorBrewer)
library(ggeffects)
library(ggridges)
library(ggrepel)
library(tidylog)
library(gt)
library(patchwork)
library(sf)
library(mapview)
library(basemaps)
library(ggtext)
library(rlang)
library(tidyverse)

# Load data ---------------------------------------------------------------

##IPM Output
out <- readRDS(file = here::here("jags/output/BCAB_CaribouIPM_posteriors.rds"))

##IPM input to compare results
hd <- read.csv("data/clean/blueprint.csv")
hn <- hd %>%
	dplyr::select(herd, herd_num)
afs <- read.csv("data/clean/survival.csv") %>%
	arrange(herd) %>%
	left_join(hn, by = "herd")
afr <- read.csv("data/clean/recruitment.csv") %>%
	arrange(herd) %>%
	left_join(hn, by = "herd")
counts <- read.csv("data/clean/counts.csv") %>%
	arrange(herd) %>%
	left_join(hn, by = "herd")
trt <- read.csv("data/clean/treatments.csv") %>%
	arrange(herd) %>%
	left_join(hn, by = "herd")
ecotype <- read.csv("data/raw/treatment.csv") %>%
  dplyr::select(herd=Herd,ECCC=ECCC_Recov_Grp, COSEWIC=COSEWIC_Grp, Heard_Vagt1998=Heard.and.Vagt.1998.grouping)%>%
  distinct()

#  Years of study
yrs <-  seq(from = min(trt$year), to = max(trt$year), by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
yr_df <- as.data.frame(cbind(yrs, yr_idx))


## ----Data housekeeping, results='hide', message=FALSE, warning=FALSE---------------------------------------------------------------------------------------
##set up colors for plotting
display.brewer.pal(8,"Accent")
cols <- RColorBrewer::brewer.pal(8, "Accent")

###Change Narraway BC to Bearhole Redwillow
trt <- trt%>%mutate(herd=case_when(herd%in%"Narraway BC"~"Bearhole Redwillow", TRUE~herd))
hd <- hd%>%mutate(herd=case_when(herd%in%"Narraway BC"~"Bearhole Redwillow", TRUE~herd))
counts <- counts%>%mutate(herd=case_when(herd%in%"Narraway BC"~"Bearhole Redwillow", TRUE~herd))
afr <- afr%>%mutate(herd=case_when(herd%in%"Narraway BC"~"Bearhole Redwillow", TRUE~herd))
afs <- afs%>%mutate(herd=case_when(herd%in%"Narraway BC"~"Bearhole Redwillow", TRUE~herd))
ecotype <- ecotype%>%mutate(herd=case_when(herd%in%"Narraway BC"~"Bearhole Redwillow", TRUE~herd))

# Pull out some values, summarize to demog data frame ---------------------

####TREATMENT COMBOS
treatment.combos <- trt%>%filter(applied%in%1)%>%
  select(-intensity)%>%
  pivot_wider(names_from=treatment, values_from=treatment)%>%
  mutate(trt=paste(`reduce wolves`,`sterilize wolves`,`reduce moose`,pen,feed,transplant, sep="-")%>%
           str_replace_all(" ","")%>%
           str_replace_all("-NA","")%>%
           str_replace_all("NA-",""))%>%
  select(herd,yrs=year,trt)

##pull posterior draws, add in herd, year, and treatment to each herd-year
demog <- out %>%
  spread_draws(c(totNMF, totN,totAdults, S,R_adj, lambda, R.ad, SR)[i,j]) %>%
  median_qi(.width = 0.9)%>%
  left_join(yr_df%>%rename(j=yr_idx))%>%
  left_join(hd%>%select(herd,i=herd_num))%>%
  left_join(treatment.combos)%>%
  mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))%>%
  rename(R=R_adj, R.lower=R_adj.lower, R.upper=R_adj.upper)


##filter herds to a period where not fully extirpated
extirpated.yr <- demog%>%
  filter(round(totN,0)<=1)%>%
  group_by(herd)%>%
  filter(yrs==min(yrs))%>%
  ungroup


herds<- unique(demog$herd)

demog.trim <- tibble()
for(i in 1:length(herds)){
  
  if(herds[i]%in%extirpated.yr$herd){
  yr <-extirpated.yr%>%
    dplyr::filter(herd==!!herds[i])
  a <- demog%>%
    dplyr::filter(herd==!!herds[i] & yrs<yr$yrs)}
  
  if(!herds[i]%in%extirpated.yr$herd){
    a <- demog%>%
      dplyr::filter(herd==!!herds[i])}
  
  demog.trim <- bind_rows(a,demog.trim)
}

demog <- demog.trim
rm(demog.trim)
rm(a)
rm(yr)

write_csv(demog,"tables/demog.csv")

demog.mod <- demog%>%filter(totAdults>10 & totNMF>20) ##modelling data that doesn't include functionally extirpated herds, demography gets unstable


## ----Plot herd abundance, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------

##Prep data and layout for plot
trt.plot <-trt%>%
  filter(applied==1)%>%
  left_join(demog%>%group_by(herd)%>%summarize(max=max(totNMF.upper)))%>%
  mutate(y=case_when(treatment%in%"reduce wolves"~max-(max*0.1),
                     treatment%in%"transplant"~max-(max*0.15),
                     treatment%in%"reduce moose"~max-(max*0.2),
                     treatment%in%"sterilize wolves"~max-(max*0.25),
                     treatment%in%"pen"~max-(max*0.3),
                     treatment%in%"feed"~max-(max*0.35)))

element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
    element$family <- element$hi.family %||% element$family
  }
  NextMethod()
}

sims <- out %>%
  gather_draws(c(pred_totNMF,totNMF,totAdults, totCalvesMF,pred_totCalvesMF)[i,j]) %>%
  median_qi(.width = 0.9)%>%
  left_join(yr_df%>%rename(j=yr_idx))%>%
  left_join(hd%>%select(herd,i=herd_num))%>%
  left_join(treatment.combos)%>%
  mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))

sims.plot <- sims%>%
  filter(!.variable%in%c("totAdults","totCalvesMF","pred_totCalvesMF"))%>%
  group_by(yrs,.variable)%>%
  summarise(across(.value:.upper,~sum(.x)))%>%
  mutate(.variable=case_when(.variable%in%"pred_totNMF"~"Status quo",
                             TRUE~"Observed"))%>%
  mutate(.variable=fct_relevel(.variable,"Status quo","Observed"))

ext.yr <- sims%>%
  filter(.variable%in%c("totAdults","totNMF"))%>%
  select(.variable,.value,yrs,herd)%>%
  pivot_wider(names_from=".variable", values_from=".value")%>%
  group_by(herd,yrs)%>%
  summarise(totAdults=median(totAdults),
            totNMF=median(totNMF))%>%
  group_by(herd)%>%
  mutate(viable=case_when(yrs==2021 & round(totAdults,0)>8 & round(totNMF,0)>20 ~1, TRUE~0),
         viable=max(viable))%>%
  filter(round(totAdults,0)<=8 | round(totNMF,0)<=20, viable==0)%>% ##pull if 8 or fewer females present, or 20 or fewer total
  filter(yrs==(min(yrs)))%>%
  left_join(sims.plot%>%filter(.variable%in%"Observed")%>%select(yrs,.value))

##Herd Abundance
ggplot() +
  geom_ribbon(data=demog,aes(x=yrs, y=totNMF, ymin=totNMF.lower, ymax=totNMF.upper),alpha=0.4, color=NA, fill=cols[3])+
  geom_line(data=demog,aes(x=yrs, y=totNMF, ymin=totNMF.lower, ymax=totNMF.upper),size=1, color=cols[3]) +
  geom_rug(data=rbind(afr%>%select(herd,year,est)%>%mutate(type="Recruit"),afs%>%select(herd,year,est)%>%mutate(type="Surv"),counts%>%select(herd,year,est=Est_CL)%>%mutate(type="Count")),
           aes(x=year), sides = "t",length = unit(0.05, "npc"),alpha=0.5)+
  theme_ipsum()+
  theme(legend.position = "none")+
  ylab("")+
  xlab("Year")+
  facet_wrap(vars(herd), scales="free_y", ncol=8)+
  labs(x="Year",title="Abundance")+
  expand_limits(y=0)+
  scale_x_continuous(    limits = c(1974, 2026),
                         breaks = seq(1980, 2020, by = 20))+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text = element_textbox_highlight(
          size = 14,
          fill = "white", box.color = "white", 
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(4, 0, 4, 0), margin = margin(4, 1, 4, 1),
          hi.labels = ext.yr$herd, 
          hi.fill = "firebrick", hi.box.col = "firebrick", hi.col = "white"
        ),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
  geom_point(data=counts,aes(x=year, y=Est_CL),size=0.5, alpha=0.7)+
  geom_linerange(data=counts%>%mutate(Est_CL.max=case_when(Est_CL.max>5000~5000,TRUE~Est_CL.max)),aes(x=year, ymin=Est_CL.min,  ymax=Est_CL.max),alpha=0.5)+
  geom_point(data=trt.plot,aes(x=year, y=y, group=treatment, color=treatment),size=0.5)+
  scale_color_manual(values=cols[-4])+
  geom_text(data = trt.plot%>%distinct(herd,treatment,y)%>%mutate(t=str_remove(treatment,"reduce ")%>%str_sub(1,1)), aes(label = t, colour = treatment, x = 2025, y = y),
                           direction = "y")+
  coord_cartesian(clip = "off")


ggsave(here::here("plots","abundance.png"), width=15, height=11, bg="white")


## ----Plot total abundance, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------

####Summarize bou pop in '73 vs 2021####
sims.summary <- sims%>%
  group_by(herd)%>%
  filter(yrs==min(yrs)|yrs==max(yrs))%>%
  arrange(yrs)%>%
  mutate(time=case_when(yrs==min(yrs)~"1973",
                        TRUE~"2021"))%>%
  group_by(time,.variable)%>%
  summarise(abund=sum(.value),
            lower=sum(.lower),
            upper=sum(.upper))%>%
  ungroup

##what was the decline?
p.decline <- sims.summary%>%
  filter(.variable=="totNMF")%>%
  select(time, abund)%>%
  pivot_wider(names_from=time, values_from=abund)%>%
  summarise(dif=((`1973`-`2021`)/`1973`)*100)%>%
  round(0)

##how many more caribou now than status quo w/ error
sims.draws <- out %>%
  gather_draws(pred_totNMF[i,j],totNMF[i,j],totCalvesMF[i,j],pred_totCalvesMF[i,j])

n.recovery.all <-sims.draws%>%
  filter(j%in%c(49),.variable%in%c("totNMF","pred_totNMF"))%>%
  group_by(.draw,.variable)%>%
  summarise(across(.value,~sum(.x)))%>%
  pivot_wider(names_from=.variable, values_from=.value)%>%
  mutate(dif=totNMF-pred_totNMF)%>%
  pull(dif)

quantile(n.recovery.all,c(0.05,0.5,0.95))%>%round(0)

n.recovery <- median(n.recovery.all)%>%round(0)

##do again but just for calves
calves.recovered.all <-sims.draws%>%
  filter(.variable%in%c("totCalvesMF","pred_totCalvesMF"))%>%
  group_by(.draw,.variable,j)%>%
  summarise(across(.value,~sum(.x)))%>%
  pivot_wider(names_from=.variable, values_from=.value)%>%
  mutate(dif=totCalvesMF-pred_totCalvesMF)%>%
  group_by(.draw)%>%
  summarise(dif=sum(dif))%>%
  pull(dif)

quantile(calves.recovered.all,c(0.05,0.5,0.95))%>%round(0)

calves.recovered <- median(calves.recovered.all)%>%round(0)
  
sum(out$mean$totCalvesMF-out$mean$pred_totCalvesMF)


##with error


median(n.recovery.all)
quantile(n.recovery.all,c(0.05,0.5,0.95))


####Total Abundance####
abundance.all.plot <- ggplot(data=sims.plot) +
  #ggplot(aes(x=yrs,  y=.value, color=.variable,fill=.variable)) +
    geom_ribbon(alpha=0.3, aes(x=yrs,  y=.value, ymin=.lower, ymax=.upper, color=.variable,fill=.variable))+
    geom_line(size=1, aes(x=yrs,  y=.value, ymin=.lower, ymax=.upper, color=.variable,fill=.variable)) +
    geom_text(data = sims.plot%>%filter(yrs==2021), aes(label = .variable, colour = .variable, x = Inf, y = .value), hjust = 0) +
    geom_jitter(data=ext.yr,size=1, aes(x=yrs,  y=.value),alpha=0.5) +
    theme_ipsum()+
    theme(legend.position = "none")+
    ylab("")+
    xlab("Year")+
    labs(x="Year", y="Abundance",title="b) Population Trend",
         subtitle=paste0(p.decline,"% decline since 1973, +", n.recovery," caribou from recovery measures"))+
    expand_limits(y=0)+
   scale_y_continuous(expand = c(0, 100))+
    scale_x_continuous(breaks = seq(1980, 2020, by = 20))+
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          strip.text.x = element_text(size=15),
          strip.text.y = element_text(size=15),
          legend.text = element_text(size=13),
          legend.title=element_text(size=15),
          plot.margin = unit(c(1,5,1,1), "lines"))+
   geom_rug(data=trt%>%filter(applied%in%1), aes(x=year), sides = "b",length = unit(0.05, "npc"),alpha=0.05)+
   scale_color_manual(values=cols[c(3,1)])+
  annotate(geom = "text", x = 1974, y = 2200, label = "Recovery\nmeasures", hjust = "left")+
  annotate(geom = "curve", x = 1976, y = 1200, xend = 1985, yend = 300, 
    curvature = .3, arrow = arrow(length = unit(2, "mm")))+
  annotate(geom = "text", x = 1985, y = 8500, label = "Extirpation\nevent", hjust = "left")+
  annotate(geom = "curve", x = 1985, y = 8500, xend = (ext.yr%>%ungroup%>%filter(yrs==min(yrs))%>%pull(yrs))+1, yend = (ext.yr%>%ungroup%>%filter(yrs==min(yrs))%>%pull(.value))+100, 
           curvature = .3, arrow = arrow(length = unit(2, "mm")))+
  coord_cartesian(
    clip = "off"
  )
abundance.all.plot

ggsave(plot=abundance.all.plot,here::here("plots","abundance_all.png"), width=6, height=6, bg="white")


## ----Plot total abundance ecotype, message=FALSE, warning=FALSE,eval=FALSE,include=FALSE-------------------------------------------------------------------
## ##for AB only:filter(!herd%in%c("A La Peche","Banff","Brazeau","Maligne","Narraway AB","Tonquin","Redrock/Prairie Creek"))%>%
## 
## sims.plot.grp <- sims%>%
##   filter(.variable%in%c("pred_totNMF","totNMF"))%>%
##   left_join(ecotype)%>%
##   group_by(yrs,.variable,ECCC)%>%
##   summarise(across(.value:.upper,~sum(.x)))%>%
##   mutate(.variable=case_when(.variable%in%"pred_totNMF"~"Status quo*",
##                              TRUE~"Observed"))%>%
##   ungroup%>%
##   mutate(.variable=fct_relevel(.variable,"Status quo*","Observed"),
##          ECCC=fct_relevel(ECCC,"Southern Group","Central Group", "Northern Group"))
## 
## ext.yr.grp <- sims%>%
##   filter(.variable=="totNMF")%>%
##   group_by(herd,yrs)%>%
##   summarise(abund=median(.value))%>%
##   filter(abund<20)%>%
##   filter(yrs==(min(yrs)))%>%
##   left_join(ecotype)%>%
##   left_join(sims.plot.grp%>%filter(.variable%in%"Observed")%>%select(yrs,.value,ECCC))%>%
##   select(herd,yrs,abund,ECCC,.value)%>%
##   ungroup%>%
##   mutate(ECCC=fct_relevel(ECCC,"Southern Group","Central Group", "Northern Group"))
## 
## 
## ggplot(data=sims.plot.grp)+
##   geom_ribbon(alpha=0.3, aes(x=yrs,  y=.value, ymin=.lower, ymax=.upper, color=.variable,fill=.variable))+
##   geom_path(size=1, aes(x=yrs,  y=.value, ymin=.lower, ymax=.upper, color=.variable,fill=.variable)) +
##   geom_text(data = sims.plot.grp%>%filter(yrs==2021), aes(label = .variable, colour = .variable, x = Inf, y = .value), hjust = 0) +  theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   labs(x="Year", y="Abundance",title="SMC Population Trend",
##        subtitle=paste0("",p.decline,"% decline since 1974, +", n.recovery," caribou from recovery measures"),
##        caption ="*Predicted population trend under a status quo scenario where no recovery measures were implemented")+
##   expand_limits(y=0)+
##   scale_y_continuous(expand = c(0, 100))+
##   scale_x_continuous(breaks = seq(1980, 2020, by = 20))+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=12),
##         strip.text.y = element_text(size=15),
##         plot.caption = element_text(color=cols[c(3)],size=12),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15),
##         plot.margin = unit(c(1,5,1,1), "lines"))+
##   geom_rug(data=trt%>%left_join(ecotype)%>%filter(applied%in%1), aes(x=year), sides = "b",length = unit(0.05, "npc"),alpha=0.05)+
##   scale_color_manual(values=cols[c(3,1)])+
##   geom_text(data=tibble(lab="Recovery\nmeasures", yrs = 1990, .value = 1300,ECCC="Northern Group"), aes(x=yrs,  y=.value, label="Recovery\nmeasures"),hjust = "left")+
##   geom_curve(data=tibble(x = 2002, y = 1000, xend = 2010, yend = 500,ECCC="Northern Group"),
##              aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE,curvature = -.3, arrow = arrow(length = unit(2, "mm")))+
##   coord_cartesian(
##     clip = "off"
##   )+
##   facet_wrap(vars(ECCC),scales="free_y")
## 
## 
## ggsave(here::here("plots","abundance_bygroup.png"), width=13, height=6, bg="white")


## ----Plot vital rates, message=FALSE, warning=FALSE,eval=FALSE,include=FALSE-------------------------------------------------------------------------------
## ####Survival####
## ggplot() +
##   geom_ribbon(data=demog,aes(x=yrs, y=S, ymin=S.lower, ymax=S.upper, fill=herd),alpha=0.4)+
##   geom_point(data=demog,aes(x=yrs, y=S, ymin=S.lower, ymax=S.upper, fill=herd),alpha=0.7) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Survival")+
##   #expand_limits(y=0)+
##   scale_x_continuous(breaks = seq(1980, 2020, by = 20))+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_point(data=afs,aes(x=year, y=est, color=herd),alpha=0.7)+
##   geom_linerange(data=afs%>%mutate(ucl=case_when((est+sd)>1~1,
##                                                  TRUE~est+sd),
##                                    lcl=case_when((est-sd)<0~0,
##                                                  TRUE~est-sd)),aes(x=year, ymin=lcl,  ymax=ucl, color=herd),alpha=0.7)
## ggsave(here::here("plots","survival.png"), width=14, height=11, bg="white")
## 
## ####Recruitment####
## ggplot() +
##   geom_ribbon(data=demog,aes(x=yrs, y=R, ymin=R.lower, ymax=R.upper, fill=herd),alpha=0.4)+
##   geom_line(data=demog,aes(x=yrs, y=R, ymin=R.lower, ymax=R.upper, fill=herd)) +
##   geom_point(data=demog,aes(x=yrs, y=R, ymin=R.lower, ymax=R.upper, fill=herd)) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   #expand_limits(y=0)+
##   scale_x_continuous(breaks = seq(1980, 2020, by = 20))+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Recruitment")+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_point(data=afr,aes(x=year, y=est, color=herd))+
##   geom_linerange(data=afr%>%mutate(ucl=case_when((est+sd)>1~1,
##                                                  TRUE~est+sd),
##                                    lcl=case_when((est-sd)<0~0,
##                                                  TRUE~est-sd)),aes(x=year, ymin=lcl,  ymax=ucl, color=herd),alpha=0.7)
## ggsave(here::here("plots","recruitment.png"), width=14, height=11, bg="white")
## 
## 
## #write_csv(demog.trim, here::here("tables","demog.trim.csv"))
## 
## 
## ####Adult Recruitment####
## ggplot() +
## geom_ribbon(data=demog,aes(x=yrs, y=R.ad, ymin=R.ad.lower, ymax=R.ad.upper, fill=herd),alpha=0.4)+
##   geom_line(data=demog,aes(x=yrs, y=R.ad, ymin=R.ad.lower, ymax=R.ad.upper, fill=herd)) +
##   geom_point(data=demog,aes(x=yrs, y=R.ad, ymin=R.ad.lower, ymax=R.ad.upper, fill=herd)) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Recruitment (calves/adult)")+
##   expand_limits(y=0)+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))
## 
## ####Lambda####
## ggplot() +
##   geom_ribbon(data=demog,aes(x=yrs, y=lambda, ymin=lambda.lower, ymax=lambda.upper, fill=herd),alpha=0.4)+
##   geom_line(data=demog,aes(x=yrs, y=lambda, ymin=lambda.lower, ymax=lambda.upper, fill=herd)) +
##   geom_point(data=demog,aes(x=yrs, y=lambda, ymin=lambda.lower, ymax=lambda.upper, fill=herd)) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Population Growth")+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_hline(data=data.frame(herd=unique(demog$herd)),
##              aes(yintercept =1),linetype="dashed")
## 


## ----trt eff- r, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------

##Gather draws
demog.draws <- out %>%
  gather_draws(lambda[i,j],S[i,j],R_adj[i,j],R.ad[i,j],totNMF[i,j])%>%
  mutate(.variable=case_when(.variable%in%"R_adj"~"R",
                             TRUE~.variable))

###pivot wider
demog.draws  <- demog.draws%>%
  ungroup%>%
  pivot_wider(names_from=.variable, values_from=.value)
  

##add in years and herd names
demog.draws  <- demog.draws %>%
  left_join(yr_df%>%rename(j=yr_idx), by="j")%>%
  left_join(trt%>%distinct(herd_num, herd)%>%select(i=herd_num,herd),by=c("i"))%>%
  left_join(treatment.combos)%>%
  left_join(trt%>%distinct(herd, year, intensity)%>%select(herd,yrs=year, intensity)%>%filter(intensity%in%"low"),by=c("herd","yrs"))%>%
  mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))


##remove first year lambda for each herd, as lambda==1
demog.draws  <- demog.draws %>%
  group_by(herd, .draw)%>%
  arrange(yrs)%>%
  slice(-1)%>%
  ungroup


##filter herds to a period where not extirpated
demog.draws.trim <- tibble()
for(i in 1:length(herds)){
  
  if(herds[i]%in%ext.yr$herd){
    yr <-ext.yr%>%
      dplyr::filter(herd==!!herds[i])
    a <- demog.draws %>%
      dplyr::filter(herd==!!herds[i] & yrs<yr$yrs)}
  
  if(!herds[i]%in%ext.yr$herd){
    a <- demog.draws %>%
      dplyr::filter(herd==!!herds[i])}
  
  demog.draws.trim <- bind_rows(a,demog.draws.trim)
}

demog.draws  <- demog.draws.trim
rm(demog.draws.trim)


###Remove years where transplant was being done (impossibly high lambda due to adding individuals, not population response)
demog.draws  <- demog.draws %>%
  filter(!(herd%in%"Charlotte Alplands" & yrs%in%c(1984:1991)),
         !(herd%in%"Telkwa" & yrs%in%c(1997:1999)),
         !(herd%in%"South Selkirks" & yrs%in%c(1987:1989)),
         !(herd%in%"Purcell South" & yrs%in%c(2012:2013)))

# ###how to calculate from draws
# lambda.compare <- demog%>%
#   ungroup%>%
#   dplyr::select(herd,yrs,lambda)%>%
#   left_join(demog.draws %>%
#               group_by(herd,yrs)%>%
#               summarise(median=median(lambda),
#                         mean=mean(lambda)))
# 
# lm(lambda~median, data=lambda.compare)%>%summary
# lm(lambda~mean, data=lambda.compare)%>%summary
# 
# ##definitely median when calculating for each posterior. Geo mean likely more appropriate for summarizing pooled distributions (i.e., overall means)

demog.draws.combotreat <- demog.draws %>%
  #filter(!(intensity%in%"low"|totNMF<20))%>%
  group_by(.draw,trt,herd)%>%
  summarise(lambda=exp(mean(log(lambda))))%>% ##geo mean per herd-treatment-draw
  group_by(.draw,trt)%>%
  summarise(lambda=exp(mean(log(lambda))))%>% ##geo mean lambda per treatment-draw
  mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))%>%
  mutate(group=case_when(trt%in%"Reference"~"Reference", TRUE~"Treatment"))

demog.draws.combotreat.rug <- demog.draws %>%
  #filter(!(intensity%in%"low"|totNMF<20))%>%
  group_by(trt,herd)%>%
  summarise(lambda=exp(mean(log(lambda))))%>% ##geo mean per herd-treatment-draw
  mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))%>%
  mutate(group=case_when(trt%in%"Reference"~"Reference", TRUE~"Treatment"))

order <- demog.draws.combotreat%>%
  group_by(trt)%>%
  summarise(med=median(lambda))

demog.draws.combotreat%>%
  left_join(order)%>%
  filter(trt!="transplant")%>%
    ggplot(aes(x=log(lambda), y=fct_reorder(trt,med), fill=group))+
  geom_density_ridges(#scale = 1.5, 
    #scale = 1.3, 
    rel_min_height = .01,
    size = 0.25,
    alpha=0.9)+
    theme_ipsum()+
  geom_point(data=demog.draws.combotreat.rug%>%
               left_join(order)%>%
               filter(trt!="transplant"),
             aes(y=fct_reorder(trt,med),x=log(lambda)),
             shape = "|")+
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          strip.text.x = element_text(size=15),
          strip.text.y = element_text(size=15),
          legend.text = element_text(size=13),
          legend.title=element_text(size=15),
          legend.position = "none")+
    geom_vline(xintercept = 0, linetype="dashed")+
   xlim(-0.25,0.3)+
  labs(x="Instantaneous rate of increase (r)",
       y="Recovery measure(s)",
       fill="",
       title="Instantaneous Rate of Increase by Recovery Measure")+
  scale_fill_manual(values=cols[c(3,1)])

ggsave(here::here("plots","lambda_treatments.png"), width=8, height=7, bg="white")




demog.draws.combotreat%>%
  group_by(trt)%>%
  summarise(lambda.med=median(lambda,na.rm=TRUE)%>%round(2),
            lower=quantile(lambda,0.05,na.rm=TRUE)%>%round(2),
            upper=quantile(lambda,0.95,na.rm=TRUE)%>%round(2))%>%
  arrange(-lambda.med)%>%
  mutate(lambda=paste0(lambda.med," (",lower,"-",upper,")"))%>%
  gt()


## ----trt eff- BA, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------------
# Before-After ---------------------------------------------

##pull draws and organize into treatment and untreated (reference) timeframes for each herd
eff.draws <- demog.draws%>%
  mutate(lambda=log(lambda))%>% ##little r
  ungroup()%>%
  filter(!trt%in%"Reference")%>%
  group_by(.draw,herd,trt)%>%
  summarise(across(lambda:R, ~mean(.x)))%>% 
  pivot_longer(lambda:R)%>%
  ungroup()%>%
  select(.draw,trt,herd,name,eff=value)%>%
  left_join(demog.draws%>%
              mutate(lambda=log(lambda))%>% ##little r
              ungroup()%>%
              filter(trt%in%"Reference")%>%
              group_by(.draw,herd)%>%
              summarise(across(lambda:R,  ~mean(.x)))%>% 
              pivot_longer(lambda:R)%>%
              select(.draw,herd,name,ref=value),
            by=c(".draw","herd","name"))%>%
  mutate(delta.i=exp(eff-ref)-1,
         delta.i.lr=eff-ref)



order <- eff.draws%>%
  filter(name=="lambda")%>%
  group_by(trt)%>%
  summarise(med=median(delta.i.lr))%>%
  arrange(-med)

eff.draws <-eff.draws%>%
  left_join(order)%>%
  mutate(trt=fct_reorder(trt,med),
         name=case_when(name=="lambda"~"Rate of increase (r)",
                        name=="R"~"Recruitment (R)",
                        name=="S"~"Survival (S)",
                        TRUE~name))


ggplot()+
  geom_density_ridges(data= eff.draws%>%
                        filter(trt!="transplant")%>%
                        group_by(name,trt,.draw)%>%
                        summarise(delta=median(delta.i.lr)),
                      aes(x=delta, y=trt,fill=name),
                      scale = .9, 
                      rel_min_height = .01,
                      size = 0.25,
                      alpha=0.9)+
  geom_point(data=eff.draws%>%
               filter(trt!="transplant")%>%
               group_by(herd,name,trt)%>%
               summarise(delta=median(delta.i.lr)),
             aes(y=trt,x=delta),
             shape = "|")+
  theme_ipsum()+
  facet_wrap(vars(name))+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position="none")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="Change in value", y="Recovery measure(s)", title="Before-After Assessment of Effectivness")+
  xlim(-0.2,0.4)+
  scale_fill_manual(values=cols[c(1:3)])

ggsave(here::here("plots","ba_all.png"), width=10, height=7, bg="white")


trt_eff_ba_table <- eff.draws%>%
  filter(name=="Rate of increase (r)")%>%
  group_by(trt,.draw)%>%
  summarise(delta=median(delta.i))%>%
  summarise(delta.l=median(delta,na.rm=TRUE)%>%round(2),
            lower=quantile(delta,0.05,na.rm=TRUE)%>%round(2),
            upper=quantile(delta,0.95,na.rm=TRUE)%>%round(2))%>%
  arrange(-delta.l)%>%
  mutate(delta.lambda=paste0(delta.l," (",lower,"-",upper,")"))%>%
  select(Treatment=trt,delta.lambda)%>%
  filter(Treatment!="transplant")

trt_eff_ba_table%>%
  gt()%>%
  gtsave(here::here("tables","trt_eff_ba.rtf"))


eff.draws%>%
  filter(name=="Rate of increase (r)")%>%
  group_by(herd,name,trt)%>%
  summarise(delta=median(delta.i.lr))%>%
  ungroup%>%
  count(delta>0)


## ----trt eff- BA w intensity, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------

eff.draws.app <- demog.draws%>%
  group_by(herd,yrs,trt)%>%
  mutate(totNMF.median=median(totNMF))%>% #get average pop size so popsize threshold doesnt split low/standard in some years due to draws being above/below threshold
  ungroup%>%
  mutate(lambda=log(lambda),##little r
         application=case_when(intensity=="low" | totNMF.median<30~"low", TRUE~"standard"))%>% 
  ungroup()%>%
  filter(!trt%in%"Reference")%>%
  group_by(.draw,herd,trt,application)%>%
  summarise(across(lambda:R, ~mean(.x)))%>% 
  pivot_longer(lambda:R)%>%
  ungroup()%>%
  select(.draw,application,trt,herd,name,eff=value)%>%
  left_join(demog.draws%>%
              mutate(lambda=log(lambda))%>% ##little r
              ungroup()%>%
              filter(trt%in%"Reference")%>%
              group_by(.draw,herd)%>%
              summarise(across(lambda:R,  ~mean(.x)))%>% 
              pivot_longer(lambda:R)%>%
              select(.draw,herd,name,ref=value),
            by=c(".draw","herd","name"))%>%
  mutate(delta.i=exp(eff-ref)-1,
         delta.i.lr=eff-ref)

eff.draws.app.bind <- eff.draws.app%>%
  mutate(group="all")%>%
  rbind(eff.draws.app%>%
            filter(application=="standard")%>%
            mutate(group="standard"))



order <- eff.draws.app.bind%>%
  filter(name=="lambda")%>%
  group_by(trt)%>%
  summarise(med=median(delta.i.lr))%>%
  arrange(-med)

eff.draws.app.bind <-eff.draws.app.bind%>%
  left_join(order)%>%
  mutate(trt=fct_reorder(trt,med),
         name=case_when(name=="lambda"~"Rate of increase (r)",
                        name=="R"~"Recruitment (R)",
                        name=="S"~"Survival (S)",
                        TRUE~name))


ggplot()+
  geom_density_ridges(data= eff.draws.app.bind%>%
                        filter(trt!="transplant")%>%
                        group_by(name,trt,.draw,group)%>%
                        summarise(delta=median(delta.i.lr)),
                      aes(x=delta, y=trt,fill=name,linetype=group),
                      scale = .9, 
                      rel_min_height = .01,
                      size = 0.25,
                      alpha=0.7)+
  geom_point(data=eff.draws.app.bind%>%
               filter(trt!="transplant", group=="all",!(herd=="Frisby-Boulder" & application=="standard"))%>%
               group_by(name,herd,trt,application)%>%
               summarise(delta=median(delta.i.lr)),
             aes(y=trt,x=delta, color=application),
             shape = "|", size=4)+
  theme_ipsum()+
  facet_wrap(vars(name))+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="Change in value", y="Recovery measure(s)", title="Before-After Assessment of Effectivness")+
  xlim(-0.2,0.4)+
  scale_fill_manual(values=cols[c(1:3)])+
  scale_color_manual(values=cols[c(8,6)])+
  guides(fill = "none")

ggsave(here::here("plots","ba_all_intensity.png"), width=10, height=7, bg="white")


# ggplot()+
#   geom_density_ridges(data= eff.draws.app.bind%>%
#                         filter(trt=="reducemoose",name=="Rate of increase (r)", application=="low"),
#                       aes(x=delta.i.lr, y=trt,fill=herd,linetype=herd),
#                       scale = .9, 
#                       rel_min_height = .01,
#                       size = 0.25,
#                       alpha=0.7)



eff.draws.app.bind%>%
  filter(trt!="transplant", name=="Rate of increase (r)")%>%
  group_by(.draw,name,trt,group)%>%
  summarise(delta=median(delta.i.lr,na.rm=TRUE))%>%
  group_by(trt,group)%>%
  summarise(delta.l=median(delta,na.rm=TRUE)%>%round(2),
            lower=quantile(delta,0.05,na.rm=TRUE)%>%round(2),
            upper=quantile(delta,0.95,na.rm=TRUE)%>%round(2))%>%
  arrange(-delta.l)%>%
  arrange(trt,group)


## ----trt eff- BACI w ecotype, message=FALSE, warning=FALSE,eval=FALSE,include=FALSE------------------------------------------------------------------------
## ####ECCC grouping
## demog.mod.baci <- demog.draws%>%
##   ungroup%>%
##   left_join(ecotype)%>%
##   select(.draw,herd,yrs,lambda,totNMF,trt,ECCC)
## 
## 
## trt.herds <- demog.mod.baci%>%
##   filter(!trt%in%c("Reference","transplant"),.draw==1)%>%
##   group_by(herd,trt,ECCC)%>%
##   summarise(trt.start=min(yrs),
##             trt.end=max(yrs),
##             dur=n(),
##             baci.start=trt.start-dur)%>%
##   ungroup
## 
## baci <- tibble()
## for(i in 1:nrow(trt.herds)){
## 
##   post.trt <- demog.mod.baci%>%
##     dplyr::filter(herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c(trt.herds[i,"trt.start"][[1]]:trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=1)
## 
##   pre.trt <- demog.mod.baci%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):(trt.herds[i,"trt.start"][[1]]-1)))%>%
##     dplyr::mutate(CI=1)
## 
##   control <-demog.mod.baci%>%
##     ungroup%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   !herd%in%!!trt.herds[i,"herd"][[1]],
##                   ECCC%in%!!trt.herds[i,"ECCC"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=0)
## 
##   ##only keep herds with full data
## if(all(nrow(pre.trt)>0,nrow(post.trt)>0,nrow(control)>0)){
##   a <- rbind(pre.trt, post.trt,control)%>%
##     dplyr::mutate(herd.grp=trt.herds[i,"herd"][[1]],
##                   trt.grp=trt.herds[i,"trt"][[1]],
##                   BA=case_when(yrs<trt.herds[i,"trt.start"][[1]]~0,
##                                yrs>=trt.herds[i,"trt.start"][[1]]~1),
##                   BA.break=trt.herds[i,"trt.start"][[1]])
## }
## 
##   baci <- rbind(baci,a)
## }
## 
## baci <- baci%>%mutate(lambda=log(lambda))%>% ##little r
##   mutate(grp=paste(herd.grp,trt, sep="_"),
##          facet=paste(herd.grp,trt.grp, sep="_"))
## 
## ggplot() +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF))%>%
##               filter(CI==1),aes(x=yrs, y=lambda,  group=herd), color="red",alpha=1, size=1) +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF)),aes(x=yrs, y=lambda, group=herd), alpha=0.3) +
##   theme_ipsum()+
##   ylab("")+
##   xlab("Year")+
##   labs(x="Year",y="Instantaneous rate of increase (r)", title="BACI set up", subtitle="control herds shown in grey")+
##   expand_limits(y=0)+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_vline(data=baci%>%
##                distinct(BA.break, facet,herd.grp,trt.grp),aes(xintercept = BA.break-0.5), linetype="dashed") +
##   facet_wrap(vars(trt.grp,herd.grp))
## ggsave(here::here("plots", "appendix","baci_setup_ECCC.png"), width=13, height=9, bg="white")
## 
## baci.draws <- baci%>%
##   group_by(.draw,trt.grp) %>%
##   do(tidy(lmer(lambda~BA+CI+BA*CI + (1|grp), data=.)))%>%
##   filter(term%in%"BA:CI")
## 
## ggplot(data=baci.draws, aes(x=estimate, y=fct_reorder(trt.grp,estimate)))+
##   geom_density_ridges(scale = .9,
##                       rel_min_height = .01,
##                       size = 0.25,
##                       alpha=0.9,
##                       fill=cols[c(1)])+
##   theme_ipsum()+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15),
##         legend.position="none")+
##   geom_vline(xintercept = 0, linetype="dashed")+
##   labs(x="Change in rate of increase", y="Recovery measure(s)", title="Before-After-Control-Impact", subtitle="w/ spatio-temporal matching using ECCC ecotype")+
##   xlim(-0.3,0.6)
## ggsave(here::here("plots", "appendix","baci_all_ECCC.png"), width=10, height=7, bg="white")
## 
## baci.draws%>%
##   group_by(trt.grp)%>%
##   summarise(delta.l=median(estimate,na.rm=TRUE)%>%round(2),
##             lower=quantile(estimate,0.05,na.rm=TRUE)%>%round(2),
##             upper=quantile(estimate,0.95,na.rm=TRUE)%>%round(2))%>%
##   arrange(-delta.l)%>%
##   mutate(delta.lambda=paste0(delta.l," (",lower,"-",upper,")"))%>%
##   select(Treatment=trt.grp,`Delta r (BACI)`=delta.lambda)%>%
##   left_join(trt_eff_ba_table%>% select(Treatment,`Delta r (BA)`=delta.lambda))%>%
##   gt()%>%
##  gtsave(here::here("tables", "appendix","trt_eff_baci_ECCC.rtf"))
## 
## 
## 
## 
## ####Heard and Vagt grouping
## demog.mod.baci <- demog.draws%>%
##   ungroup%>%
##   left_join(ecotype)%>%
##   select(.draw,herd,yrs,lambda,totNMF,trt,Heard_Vagt1998)
## 
## 
## trt.herds <- demog.mod.baci%>%
##   filter(!trt%in%c("Reference","transplant"),.draw==1)%>%
##   group_by(herd,trt,Heard_Vagt1998)%>%
##   summarise(trt.start=min(yrs),
##             trt.end=max(yrs),
##             dur=n(),
##             baci.start=trt.start-dur)%>%
##   ungroup
## 
## baci <- tibble()
## for(i in 1:nrow(trt.herds)){
## 
##   post.trt <- demog.mod.baci%>%
##     dplyr::filter(herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c(trt.herds[i,"trt.start"][[1]]:trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=1)
## 
##   pre.trt <- demog.mod.baci%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):(trt.herds[i,"trt.start"][[1]]-1)))%>%
##     dplyr::mutate(CI=1)
## 
##   control <-demog.mod.baci%>%
##     ungroup%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   !herd%in%!!trt.herds[i,"herd"][[1]],
##                   Heard_Vagt1998%in%!!trt.herds[i,"Heard_Vagt1998"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=0)
## 
##   ##only keep herds with full data
##   if(all(nrow(pre.trt)>0,nrow(post.trt)>0,nrow(control)>0)){
##     a <- rbind(pre.trt, post.trt,control)%>%
##       dplyr::mutate(herd.grp=trt.herds[i,"herd"][[1]],
##                     trt.grp=trt.herds[i,"trt"][[1]],
##                     BA=case_when(yrs<trt.herds[i,"trt.start"][[1]]~0,
##                                  yrs>=trt.herds[i,"trt.start"][[1]]~1),
##                     BA.break=trt.herds[i,"trt.start"][[1]])
##   }
## 
##   baci <- rbind(baci,a)
## }
## 
## baci <- baci%>%mutate(lambda=log(lambda))%>% ##little r
##   mutate(grp=paste(herd.grp,trt, sep="_"),
##          facet=paste(herd.grp,trt.grp, sep="_"))
## 
## ggplot() +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF))%>%
##               filter(CI==1),aes(x=yrs, y=lambda,  group=herd), color="red",alpha=1, size=1) +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF)),aes(x=yrs, y=lambda, group=herd), alpha=0.3) +
##   theme_ipsum()+
##   ylab("")+
##   xlab("Year")+
##   labs(x="Year",y="Instantaneous rate of increase (r)", title="BACI set up", subtitle="control herds shown in grey")+
##   expand_limits(y=0)+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_vline(data=baci%>%
##                distinct(BA.break, facet,herd.grp,trt.grp),aes(xintercept = BA.break-0.5), linetype="dashed") +
##   facet_wrap(vars(trt.grp,herd.grp))
## ggsave(here::here("plots", "appendix","baci_setup_Heard_Vagt1998.png"), width=13, height=9, bg="white")
## 
## baci.draws <- baci%>%
##   group_by(.draw,trt.grp) %>%
##   do(tidy(lmer(lambda~BA+CI+BA*CI + (1|grp), data=.)))%>%
##   filter(term%in%"BA:CI")
## 
## ggplot(data=baci.draws, aes(x=estimate, y=fct_reorder(trt.grp,estimate)))+
##   geom_density_ridges(scale = .9,
##                       rel_min_height = .01,
##                       size = 0.25,
##                       alpha=0.9,
##                       fill=cols[c(1)])+
##   theme_ipsum()+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15),
##         legend.position="none")+
##   geom_vline(xintercept = 0, linetype="dashed")+
##   labs(x="Change in rate of increase", y="Recovery measure(s)", title="Before-After-Control-Impact", subtitle="w/ spatio-temporal matching using Heard & Vagt 1998")+
##   xlim(-0.3,0.6)
## ggsave(here::here("plots", "appendix","baci_all_Heard_Vagt1998.png"), width=10, height=7, bg="white")
## 
## baci.draws%>%
##   group_by(trt.grp)%>%
##   summarise(delta.l=median(estimate,na.rm=TRUE)%>%round(2),
##             lower=quantile(estimate,0.05,na.rm=TRUE)%>%round(2),
##             upper=quantile(estimate,0.95,na.rm=TRUE)%>%round(2))%>%
##   arrange(-delta.l)%>%
##   mutate(delta.lambda=paste0(delta.l," (",lower,"-",upper,")"))%>%
##   select(Treatment=trt.grp,`Delta r (BACI)`=delta.lambda)%>%
##   left_join(trt_eff_ba_table%>% select(Treatment,`Delta r (BA)`=delta.lambda))%>%
##   gt()%>%
##   gtsave(here::here("tables","trt_eff_baci_Heard_Vagt1998.rtf"))


## ----trt eff- BACI w intensity, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------
##prep data with individual treatments 1/0
eff.draws <- eff.draws%>%
  mutate(reducewolves=case_when(str_detect(trt,"reducewolves")~1, TRUE~0),
         sterilizewolves=case_when(str_detect(trt,"sterilizewolves")~1, TRUE~0),
         reducemoose=case_when(str_detect(trt,"reducemoose")~1, TRUE~0),
         pen=case_when(str_detect(trt,"pen")~1, TRUE~0),
         feed=case_when(str_detect(trt,"feed")~1, TRUE~0),
         transplant=case_when(str_detect(trt,"transplant")~1, TRUE~0))

ind.eff <- eff.draws%>%
  filter(name=="Rate of increase (r)")%>%
  group_by(.draw) %>%
  do(tidy(lm(delta.i.lr~reducewolves+sterilizewolves+reducemoose+pen+feed, data=.)))%>%
  filter(!term%in%"(Intercept)")
  

ind.eff%>%group_by(term)%>%
  summarise(eff=median(estimate))


order <- ind.eff%>%
  group_by(term)%>%
  summarize(med=median(estimate))

ind.eff.plot <- ggplot(ind.eff%>%left_join(order), aes(x=estimate,y=fct_reorder(term,med),fill=term))+
  geom_density_ridges(scale = .9, 
                      rel_min_height = .01,
                      size = 0.25,
                      alpha=0.5)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position="none")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="Change in rate of increase", y="Recovery measure(s)", title="a) Individual Treatment Effects",subtitle="Partitioned using regression analysis, assuming effects are additive")+
  xlim(-0.19,0.22)+
  scale_fill_manual(values=cols[c(1:6)])

#ggsave(here::here("plots","ind_effects.png"), width=5, height=6, bg="white")


ind.eff%>%
  group_by(term)%>%
  summarise(delta.l=median(estimate,na.rm=TRUE)%>%round(2),
            lower=quantile(estimate,0.05,na.rm=TRUE)%>%round(2),
            upper=quantile(estimate,0.95,na.rm=TRUE)%>%round(2))%>%
  arrange(-delta.l)%>%
  mutate(delta.lambda=paste0(delta.l," (",lower,"-",upper,")"))%>%
  select(Treatment=term,delta.lambda)%>%
  gt()


## ----cons sims, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------
n.sims    <- 1000
start.pop <-100
#sim.ref   <- demog.draws.combotreat%>%ungroup%>%filter(trt=="Reference")%>%filter(lambda<quantile(lambda,0.95),lambda>quantile(lambda,0.05))%>%sample_n(n.sims)%>%pull(lambda)
sim.ref   <- demog.draws%>%
  #filter(trt=="Reference" & totNMF<100)%>%
  filter(trt=="Reference" & totNMF<150)%>%
  group_by(.draw,trt,herd)%>%
  summarise(lambda=exp(mean(log(lambda))))%>% ##geo mean per herd-treatment-draw
  group_by(.draw,trt)%>%
  summarise(lambda=exp(mean(log(lambda))))%>% ##geo mean lambda per treatment-draw
  ungroup%>%
  filter(lambda<quantile(lambda,0.95),lambda>quantile(lambda,0.05))%>%sample_n(n.sims)%>%pull(lambda)

sim.trt   <- ind.eff%>%
  mutate(estimate=exp(estimate)-1)%>% ##back to lambda
  group_by(term)%>%
  filter(estimate<quantile(estimate,0.95),estimate>quantile(estimate,0.05))%>%
  ungroup%>%
  select(.draw,term,estimate)%>%
  pivot_wider(names_from="term",values_from="estimate")%>%
  drop_na()%>%
  sample_n(n.sims)

sim.trt$reducewolves%>%median
sim.trt$feed%>%median

year.end <- 14
sim.df    <-list()
for(i in 1:n.sims){

cons.sim.i <-tibble(yr=1:(11+year.end))%>%
  dplyr::mutate(reference=c(start.pop*(sim.ref[i])^yr),
         reducewolves=case_when(yr%in%1:10~reference,
                        yr%in%11:30~c(reference[10]*(sim.ref[i]+sim.trt$reducewolves[i])^(yr-10))),
         feed=case_when(yr%in%1:10~reference,
                        yr%in%11:30~c(reference[10]*(sim.ref[i]+sim.trt$feed[i])^(yr-10))),
         pen=case_when(yr%in%1:10~reference,
                           yr%in%11:30~c(reference[10]*(sim.ref[i]+sim.trt$pen[i])^(yr-10))),
         reducemoose=case_when(yr%in%1:10~reference,
                           yr%in%11:30~c(reference[10]*(sim.ref[i]+sim.trt$reducemoose[i])^(yr-10))),
         sterilizewolves=case_when(yr%in%1:10~reference,
                           yr%in%11:30~c(reference[10]*(sim.ref[i]+sim.trt$sterilizewolves[i])^(yr-10))),
         `reducewolves+feed`=case_when(yr%in%1:10~reference,
                                           yr%in%11:30~c(reference[10]*(sim.ref[i]+sim.trt$feed[i]+sim.trt$reducewolves[i])^(yr-10))),
         `reducewolves+pen`=case_when(yr%in%1:10~reference,
                           yr%in%11:30~c(reference[10]*(sim.ref[i]+sim.trt$pen[i]+sim.trt$reducewolves[i])^(yr-10))),
         sim=1)

cons.sim.i$sim <- i
sim.df[[i]] <- cons.sim.i
}

sim.df <- do.call(rbind,sim.df)

nudge=0.15
sim.df.plot <-sim.df %>%
  pivot_longer(reference:`reducewolves+pen`)%>%
  group_by(name,sim)%>%
  mutate(inc=sum(last(value)>value[10]),
         ext=sum(last(value)<=20),
         abund=last(value))%>%
  group_by(name,yr)%>%
  summarise(median=median(value),
            se=sd(value),
            upper=quantile(abund,0.95),
            lower=quantile(abund,0.05),
            inc=(mean(inc)*100)%>%round(0),
            ext=(mean(ext)*100)%>%round(0))%>%
  filter(!(yr<10 & name!="reference"),
         yr<=(11+year.end))%>%
  mutate(trt=paste0(name," (",inc,", ",ext,", ", lower%>%round(0),"-",upper%>%round(0),")"),
         yr=yr-10)%>%
  mutate(yr.shift=case_when(name%in%"sterilizewolves" & yr==year.end~yr+(2.5*nudge),
                      name%in%"reducemoose" & yr==year.end~yr+(3.5*nudge),
                      name%in%"reference" & yr==year.end~yr-nudge,
                      name%in%"pen" & yr==year.end~yr,
                      name%in%"feed" & yr==year.end~yr+nudge,
                      name%in%"reducewolves" & yr==year.end~yr+(2*nudge),
                      name%in%"reducewolves+pen" & yr==year.end~yr+(3*nudge),
                      name%in%"reducewolves+feed" & yr==year.end~yr+(4*nudge),
                      TRUE~yr))

# a <-sim.df.plot%>%
#   filter(yr == last(yr))%>%pull()

recov.cols <- c(cols[1:3], "black",cols[5:8])
recov.sims.plot <- ggplot()+
  annotate("rect", xmin=-10, xmax=year.end+1,ymin = 0, ymax = 20,
           alpha = .1,fill = "black")+
  annotate("text", x=-4,y = 5, label="Functionally extirpated")+
  geom_line(data=sim.df.plot%>%filter(name!="reference"),aes(x=yr,y=median,color=trt))+
  geom_line(data=sim.df.plot%>%filter(name=="reference"),aes(x=yr,y=median,color=trt),size=1.5, linetype="dashed")+
  # geom_linerange(data=sim.df.plot%>%
  #                  group_by(name)%>%
  #                  filter(yr==max(yr)),aes(x=yr.shift,y=median,color=trt, ymin=median-se,ymax=median+se), alpha=0.7)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")+
  labs(x="Years since intervention", y="Abundance",title="b) Simulated Options to Avert Caribou Extirpation", subtitle="Labels = treatment (% samples increased, % samples extirpated, 90% end abundance interval)")+
  #geom_hline(yintercept = 10, linetype="dashed")+
  geom_text_repel(data=sim.df.plot%>%
                    filter(yr == last(yr)),
    aes(color = trt, label = trt, x=yr,y=median),
    size = 4,
    direction = "y",
    xlim = c(year.end+2, 35),
    hjust = 0,
    segment.size = .7,
    segment.alpha = .5,
    segment.linetype = "dotted",
    box.padding = .4,
    seed=999)+
  coord_cartesian(
    clip = "off"
  )+
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(-10, year.end+1), 
    breaks = seq(-10,year.end+1, by = 5)
  )+
  theme(plot.margin = margin(1,7,1,1.2, "cm"))+
  scale_color_manual(values=recov.cols)+
  geom_text(data=tibble(lab="Status quo", x = -3, y=90), aes(x=x,  y=y, label="Status quo"),hjust = "left",size=5)+
  geom_curve(data=tibble(x = -2, y = 80, xend = -3, yend = 60),
             aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE,curvature = -.3, arrow = arrow(length = unit(2, "mm")))
recov.sims.plot



##plot individual treatments together
recov.together <- ind.eff.plot+recov.sims.plot+plot_layout(widths = c(1.3,1.6))
ggsave(plot=recov.together, here::here("plots/recov.together.png"), width=13, height=6, bg="white")


## ----map, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------
##Prep Herd Bounds
herds <- st_read(here::here("data/Spatial/herds/u_bc_herds_2021_CL.shp"))%>%
  st_transform(3005)%>%
  select(herd=HERD_NAME)%>%
  mutate(herd=case_when(herd%in%"Narraway"~"Bearhole Redwillow",
                        herd%in%"Moberly"~"Klinse-Za",
                        herd%in%"Scott"~"Scott West",
                        herd%in%"Frisby Boulder"~"Frisby-Boulder",
                        herd%in%"Purcell Central"~"Purcells Central",
                        TRUE~herd))%>%
  rbind(st_read(here::here("data/Spatial/herds/Caribou_Range.shp"))%>%
          st_transform(3005)%>%
          select(herd=SUBUNIT)%>%
          mutate(herd=case_when(herd%in%"Narraway"~"Narraway AB",
                                herd%in%"Jasper"~"Brazeau",
                                herd%in%"Redrock-Prairie Creek"~"Redrock/Prairie Creek",
                                TRUE~herd)))%>%
  st_simplify(preserveTopology = FALSE, 
              dTolerance = 1000)

dist.herd <- read_csv("data/Spatial/disturbance/from_Emily/2022-01-27/herds_propdist_20220126_forCL.csv")%>%
  mutate(herd=case_when(herd%in%c("Hart South","Hart North")~"Hart Ranges",
                        herd%in%"Narraway BC"~"Bearhole Redwillow",
                        herd%in%"Central Selkirks"~"Nakusp",
                        TRUE~herd))%>%
  group_by(herd)%>%
  summarise(across(anthro_prop_dist:all_prop_dist,mean))
dist.herd <- dist.herd%>%
  rbind(dist.herd%>%filter(herd%in%"Nakusp")%>%mutate(herd="Duncan"))



herds <- dist.herd%>%
  dplyr::select(herd,human=anthro_prop_dist)%>%
  mutate(herd=case_when(herd%in%"Frisby Boulder"~"Frisby-Boulder",
                        herd%in%"Purcell Central"~"Purcells Central",
                        TRUE~herd))%>%
  left_join(demog.draws %>%
              filter(trt%in%c("Reference"))%>% 
              mutate(herd=case_when(herd%in%c("Hart South","Hart North")~"Hart Ranges",
                                    TRUE~herd))%>%
              group_by(herd)%>%
              filter(yrs%in%(max(yrs)-19):max(yrs))%>%
              summarise(lambda=exp(mean(log(lambda)))))%>% ##geo mean per herd
  left_join(ext.yr%>%mutate(ext=1)%>%ungroup%>%select(herd,ext)%>%rbind(tibble(herd="Scott West",ext=1)))%>%
  mutate(lambda=case_when(herd%in%"Scott West"~0.95,TRUE~lambda))%>% ##don't have Scott west. Extirpated
  left_join(herds)%>%
  drop_na(human)%>%
  mutate(human=round(human*100,0)%>%as.integer(),
         lambda.class=case_when(lambda<=0.95~"<0.95",
                                lambda>0.95 & lambda<=0.99~"0.95-0.99",
                                lambda>0.99 & lambda<1.01~"0.99-1.01",
                                lambda>=1.01 ~">1.01"),
         lambda.class2=case_when(lambda<=0.98~"declining",
                                lambda>0.98 & lambda<=1.02~"stable",
                                lambda>1.02 ~"increasing",
                                is.na(lambda) & ext==1~"declining"))%>%
  st_as_sf()%>%
  st_transform(4326)


# mapview(herds,zcol="human")
# mapview(herds,zcol="lambda.class2")

##Prep inset map
##Cities
cities <- st_read(here::here("data/Spatial/administrative/places.shp"))%>%
  filter(NAME%in%c("C", "VANCOUVER",  "Prince George","Fort St John"))%>%
  mutate(Name=str_to_title(NAME))%>%
  select(Name,geometry)%>%
  st_transform(4326)%>%
  rbind(tibble(Name="Banff", Y=51.18, X=-115.56)%>%st_as_sf(coords=c("X","Y"),crs=4326))%>%
  st_transform(3005)

##Pacific northwest for inset
pnw<- st_read(here::here("data/Spatial/administrative/North_America.shp"))%>%
  filter(FID_usa%in%c(1,2,8,11)|FID_canada%in%c(8,9))%>%
  mutate(prov=c("WA","MT","ID","OR","AB","BC","BC"),
         country=c(rep("USA",times=4),rep("CAN",times=3)))%>%
  st_make_valid()%>%
  group_by(prov, country)%>%
  summarise(id=mean(FID_canada))%>%
  ungroup%>%
  st_transform(3005)


library(ggspatial)
##inset
##custom crs to keep things straight
cust.crs <- "+proj=aea +lat_0=50 +lon_0=-114.9 +lat_1=49 +lat_2=50.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"
inset <- ggplot() +
  geom_sf(data = pnw%>%st_transform(cust.crs),fill="grey80", color=NA) +
  geom_sf(data = pnw%>%st_transform(cust.crs),size=0.5, fill=NA, color="grey30") +
  layer_spatial(sf::st_bbox(herds%>%st_buffer(100000)%>%st_transform(cust.crs)), fill=NA, linetype="dashed", color="grey99")+
  geom_sf_text(data=pnw%>%st_transform(cust.crs), aes(label=prov),inherit.aes = FALSE, size=3)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "transparent"),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position = c(0.65,0.075),
        plot.background = element_rect(fill = "transparent",colour = NA))




##get basemap

register_google("AIzaSyCOwGx2D77XOqRgGhKmcb5F4Kt_S61tCLI")
#set_defaults(map_service = "osm", map_type = "terrain_bg")

bmap.big <- basemap_raster(ext=herds%>%st_buffer(200000)%>%st_transform(3857),
                           map_res=1, map_type = "terrain_bg")%>%projectRaster(crs=cust.crs)


map <- ggRGB(bmap.big , r=1, g=2, b=3)+
  theme_bw()+
  geom_sf(data = pnw%>%st_transform(cust.crs),size=1, fill=NA, linetype="dashed") +
  geom_sf(data=herds,aes(fill=lambda.class2),inherit.aes = FALSE, alpha=0.7)+
  geom_sf(data=herds%>%filter(ext%in%1),aes(color="Extirpated"),fill=NA,inherit.aes = FALSE, alpha=0.7)+
  geom_sf_text(data=st_centroid(herds),aes(label = human),inherit.aes = FALSE,size=2.5, color="white")+
  geom_sf_text(data=st_centroid(herds%>%filter(lambda.class2=="stable")),aes(label = human),inherit.aes = FALSE,size=2.5, color="black")+
  geom_sf(data=cities%>%st_transform(cust.crs), inherit.aes = FALSE, size=3,pch=21, fill="white", color="black")+
  geom_sf_text(data=cities%>%filter(!Name%in%c("Fort St John","Prince George"))%>%st_transform(cust.crs), aes(label=Name),inherit.aes = FALSE, size=3,hjust = -0.1, vjust = -1, color="white")+
  geom_sf_text(data=cities%>%filter(Name%in%c("Fort St John","Prince George"))%>%st_transform(cust.crs), aes(label=Name),inherit.aes = FALSE, size=3,hjust = 0.6, vjust = -1, color="white")+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=15, color="black"),
        legend.title=element_text(size=15, color="black"),
        legend.spacing.y = unit(0.05, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.position = "top")+
  ggsn::scalebar(x.min=10E4, x.max=105E4,y.min=-19E4,y.max=85E4, dist = 150,  height=0.03, dist_unit="km",transform=FALSE, location="bottomleft", st.color = "white",st.bottom = FALSE)+
  annotation_custom(ggplotGrob(inset), xmin =65E4, xmax = 108E4, ymin = 52E4, ymax = 85E4)+
  scale_y_continuous(expand = c(0,0), limits = c(-20E4,85E4))+
  scale_x_continuous(expand = c(0,0), limits = c(5E4,105E4))+
  labs(fill="Population growth\nw/o intervention", title="a) Southern Mountain Caribou", color="")+
  scale_fill_viridis_d()+
  scale_color_manual(values = c("Extirpated"="red"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


##plot together
map_combo <- map+abundance.all.plot+plot_layout(widths = c(1.7,1.3))

ggsave(plot=map_combo, here::here("plots/map.together.png"), width=12, height=8, bg="white")




ggplot(herds%>%tibble, aes(x=human,fill=lambda.class2))+
  geom_histogram(alpha=0.5)


## ----ecotype, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------

demog.draws.combotreat.ecotype <- demog.draws %>%
  #filter(!(intensity%in%"low"|totNMF<20))%>%
  filter(yrs>=2010)%>%
  group_by(.draw,trt,herd)%>%
  summarise(lambda=exp(mean(log(lambda))))%>% ##geo mean per herd-treatment-draw
  left_join(ecotype)%>%
  group_by(.draw,trt,ECCC)%>%
  summarise(lambda=exp(mean(log(lambda))))%>% ##geo mean lambda per treatment-draw
  mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))%>%
  mutate(group=case_when(trt%in%"Reference"~"Reference", TRUE~"Treatment"))

order <- demog.draws.combotreat.ecotype %>%
  group_by(trt)%>%
  summarize(med=median(lambda))

demog.draws.combotreat.ecotype%>%
  filter(!trt%in%"transplant")%>%
  left_join(order)%>%
  ggplot(aes(x=lambda, y=fct_reorder(trt,med), fill=ECCC))+
  geom_density_ridges(#scale = 1.5, 
    #scale = 1.3, 
    rel_min_height = .01,
    size = 0.25,
    alpha=0.5)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
  geom_vline(xintercept = 1, linetype="dashed")+
  xlim(0.8,1.3)+
  labs(x="Population growth rate",
       y="Recovery measure(s)",
       fill="ECCC Ecotype",
       title="Population Growth Rate by Recovery Measure")+
  scale_fill_manual(values=cols)
ggsave(here::here("plots","appendix","lambda_treatments_ecotype.png"), width=8, height=7, bg="white")

eff.draws.ecotype <- eff.draws%>%
  dplyr::select(.draw:delta.i.lr)%>%
  filter(name=="Rate of increase (r)")%>%
  left_join(ecotype)%>%
  pivot_longer(ECCC:Heard_Vagt1998, names_to="ecotype")
  

ggplot()+
  geom_density_ridges(data= eff.draws.ecotype%>%
                        filter(!trt%in%"transplant")%>%
                        filter(ecotype=="ECCC")%>%
                        group_by(name,trt,.draw,value)%>%
                        summarise(delta=median(delta.i.lr)),
                      aes(x=delta, y=trt,fill=value),
                      scale = .9, 
                      rel_min_height = .01,
                      size = 0.25,
                      alpha=0.5)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="Delta r", y="Recovery measure(s)", title="Before-After Assessment of Effectivness", fill="ECCC Ecotype")+
  xlim(-0.2,0.4)+
  scale_fill_manual(values=cols)

ggsave(here::here("plots","appendix","baci_treatments_ecotype.png"), width=8, height=7, bg="white")

ggplot()+
  geom_density_ridges(data= eff.draws.ecotype%>%
                        filter(!trt%in%"transplant")%>%
                        filter(ecotype=="COSEWIC")%>%
                        group_by(name,trt,.draw,value)%>%
                        summarise(delta=median(delta.i.lr)),
                      aes(x=delta, y=trt,fill=value),
                      scale = .9, 
                      rel_min_height = .01,
                      size = 0.25,
                      alpha=0.5)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="delta", y="Recovery measure(s)", title="COSEWIC Grouping: Before-After Assessment of Effectivness", fill="COSEWIC Ecotype")+
  xlim(-0.2,0.4)



ggplot()+
  geom_density_ridges(data= eff.draws.ecotype%>%
                        filter(!trt%in%"transplant")%>%
                        filter(ecotype=="Heard_Vagt1998")%>%
                        group_by(name,trt,.draw,value)%>%
                        summarise(delta=median(delta.i.lr)),
                      aes(x=delta, y=trt,fill=value),
                      scale = .9, 
                      rel_min_height = .01,
                      size = 0.25,
                      alpha=0.5)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="delta", y="Recovery measure(s)", title="Before-After Assessment of Effectivness", fill="Heard_Vagt1998 Ecotype")+
  xlim(-0.2,0.4)
ggsave(here::here("plots","appendix","baci_treatments_HeardVagt1998.png"), width=8, height=7, bg="white")









ind.eff.ecotype <- eff.draws%>%
  filter(name=="Rate of increase (r)" & !trt%in%"transplant")%>%
  left_join(ecotype)%>%
  group_by(.draw) %>%
  do(tidy(lm(delta.i~reducewolves+sterilizewolves+reducemoose+pen+feed+ECCC, data=.)))%>%
  filter(!term%in%"(Intercept)")



order <- ind.eff.ecotype%>%
  group_by(term)%>%
  summarize(med=median(estimate))

ggplot(ind.eff.ecotype%>%left_join(order), aes(x=estimate,y=fct_reorder(term,med),fill=term))+
  geom_density_ridges(scale = .9, 
                      rel_min_height = .01,
                      size = 0.25,
                      alpha=0.5)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position="none")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="Delta population growth", y="Recovery measure", title="Individual Treatment Effects w/ Ecotype",subtitle="Partitioned using regression analysis, assuming effects are additive")+
  xlim(-0.15,0.25)+
  scale_fill_manual(values=cols[c(1:7)])
ggsave(here::here("plots","appendix","ind_eff_ecotype.png"), width=8, height=7, bg="white")


ind.eff.ecotype%>%group_by(term)%>%
  summarise(eff=median(estimate)%>%round(3),
            lower=quantile(estimate,0.05)%>%round(3),
            upper=quantile(estimate,0.95)%>%round(3))%>%
  arrange(-eff)%>%
  mutate(delta.lambda=paste0(eff," (",lower,"-",upper,")"))%>%
  select(Treatment=term,delta.lambda)%>%
  gt()%>%
  gtsave(here::here("tables","appendix","ind_eff_ecotype.rtf"))




ind.eff.ecotype.raw <- eff.draws%>%
  filter(name=="Rate of increase (r)" & !trt%in%"transplant")%>%
  left_join(ecotype)%>%
  dplyr::select(herd,trt,ECCC,COSEWIC,Heard_Vagt1998,reducewolves:feed,delta.i.lr)%>%
  dplyr::group_by(across(herd:feed))%>%
  summarise(delta.i.lr=median(delta.i.lr))%>%
  ungroup%>%
  mutate(ECCC=str_split(ECCC," ",simplify = TRUE)[,1])


ggplot()+
  geom_boxplot(data=ind.eff.ecotype.raw, aes(x=ECCC,y=delta.i.lr), outlier.alpha=0)+
  geom_jitter(data=ind.eff.ecotype.raw, aes(x=ECCC,y=delta.i.lr,color=trt), width=0.2)+
  geom_hline(yintercept=0,linetype="dashed")+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="ECCC Recovery Ecotype", y="Delta r")


ggplot()+
  geom_boxplot(data=ind.eff.ecotype.raw, aes(x=ECCC,y=delta.i.lr), outlier.alpha=0)+
  geom_jitter(data=ind.eff.ecotype.raw, aes(x=ECCC,y=delta.i.lr,color=trt), width=0.2)+
  geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(vars(trt))+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "none")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="ECCC Recovery Ecotype", y="Delta r")
ggsave(here::here("plots","appendix","trt_eff_boxplot_ecotype.png"), width=8, height=7, bg="white")

ggplot()+
  geom_boxplot(data=ind.eff.ecotype.raw, aes(x=COSEWIC,y=delta.i.lr), outlier.alpha=0)+
  geom_jitter(data=ind.eff.ecotype.raw, aes(x=COSEWIC,y=delta.i.lr,color=trt), width=0.2)+
  geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(vars(trt))+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "none")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="COSEWIC Recovery Ecotype", y="Delta r")

ggplot()+
  geom_boxplot(data=ind.eff.ecotype.raw, aes(x=Heard_Vagt1998,y=delta.i.lr), outlier.alpha=0)+
  geom_jitter(data=ind.eff.ecotype.raw, aes(x=Heard_Vagt1998,y=delta.i.lr,color=trt), width=0.2)+
  geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(vars(trt))+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "none")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="Heard_Vagt1998 Recovery Ecotype", y="Delta r")


ind.eff.ecotype.raw2 <- eff.draws%>%
  filter(name=="Rate of increase (r)" & !trt%in%"transplant")%>%
  left_join(ecotype)%>%
  dplyr::select(herd,trt,ECCC,COSEWIC,Heard_Vagt1998,reducewolves:feed,delta.i.lr)%>%
  dplyr::group_by(across(herd:feed))%>%
  summarise(delta.i.lr=median(delta.i.lr))%>%
  pivot_longer(reducewolves:feed)%>%
  filter(value==1)%>%
  mutate(ECCC=str_split(ECCC," ",simplify = TRUE)[,1])
ggplot()+
  geom_boxplot(data=ind.eff.ecotype.raw2, aes(x=ECCC,y=delta.i.lr), outlier.alpha=0)+
  geom_jitter(data=ind.eff.ecotype.raw2, aes(x=ECCC,y=delta.i.lr,color=trt), width=0.2)+
  geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(vars(name))+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(x="ECCC Recovery Ecotype", y="Delta r")+
  guides(colour = guide_legend(nrow = 3))
ggsave(here::here("plots","appendix","ind.trt_eff_boxplot_ecotype.png"), width=8, height=7, bg="white")


## ----for pc, message=FALSE, warning=FALSE, eval=FALSE, include=FALSE---------------------------------------------------------------------------------------
## ##For PC contract..don't review
## 
## ##Gather draws
## pc <- out %>%
##   gather_draws(S[i,j],R_adj[i,j],totAdultsMF[i,j])%>%
##   ungroup%>%
##   pivot_wider(names_from=.variable, values_from=.value)
## 
## 
## pc <- pc%>%
##   left_join(out%>%
##               gather_draws(N[i,j,a])%>%
##                ungroup%>%
##               pivot_wider(names_from=.variable, values_from=.value)%>%
##               select(i,j,.draw,a,N)%>%
##               pivot_wider(names_from=a, values_from=N,names_prefix = "N"), by=c("i","j",".draw"))%>%
##   left_join(yr_df%>%rename(j=yr_idx), by="j")%>%
##   left_join(trt%>%distinct(herd_num, herd)%>%select(i=herd_num,herd),by=c("i"))%>%
##   left_join(treatment.combos)%>%
##   left_join(trt%>%distinct(herd, year, intensity)%>%select(herd,yrs=year, intensity)%>%filter(intensity%in%"low"),by=c("herd","yrs"))%>%
##   mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))
## 
## 
## ##remove first year lambda for each herd, as lambda==1
## pc  <- pc %>%
##   group_by(herd, .draw)%>%
##   arrange(yrs)%>%
##   slice(-1)%>%
##   ungroup
## 
## 
## ##filter herds to a period where not extirpated
## pc.trim <- tibble()
## herds<- unique(demog$herd)
## for(i in 1:length(herds)){
## 
##   if(herds[i]%in%ext.yr$herd){
##     yr <-ext.yr%>%
##       dplyr::filter(herd==!!herds[i])
##     a <- pc %>%
##       dplyr::filter(herd==!!herds[i] & yrs<yr$yrs)}
## 
##   if(!herds[i]%in%ext.yr$herd){
##     a <- pc %>%
##       dplyr::filter(herd==!!herds[i])}
## 
##   pc.trim <- bind_rows(a,pc.trim)
## }
## 
## pc  <- pc.trim
## rm(pc.trim)
## 
## ##calculate N calves corrected to March
## pc <- pc%>%
##   mutate(N1_March=R_adj*totAdultsMF*0.5)
## 
## 
## ##export for Tal
## write_csv(pc%>%filter(herd%in%c("A La Peche", "Brazeau", "Maligne","Tonquin"),
##                       yrs>=2000)%>%select(.draw:intensity),"IPM_posteriors_AB.csv")
## 

