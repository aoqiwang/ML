packages <- c("shiny", "ggplot2", "reshape2", "igraph", "ggraph", "dplyr", "vegan","limma")
packages2<- "linkET"
# 检查并安装每个包
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

for (pkg in packages2) {
  if (!require(pkg, character.only = TRUE)) {
    devtools::install_github('Hy4m/linkET') 
    library(pkg, character.only = TRUE)
  }
}

library(shiny)
library(ggplot2)
library(reshape2)
library(igraph)
library(ggraph)
library(dplyr)
library(linkET)
library(vegan)

options(shiny.maxRequestSize = 500*1024^2)


ui <- fluidPage(
  titlePanel("Linkplot"),
  sidebarLayout(
    sidebarPanel(
      fileInput("datafile", "选择数据文件#行为样本列为基因(txt)", accept = c(".txt")),
      textInput("linkgene", "Link Gene", "RS"),
      textInput("heatgene", "Heatmap Genes", "GPX4,HILPDA,KIF20A,KLF2,SLC1A5,CA9"),
      numericInput("corFilter", "相关性系数阈值", value = 0),
      numericInput("pFilter", "检验P值", value = 0.05),
      numericInput("namesize", "热图字体大小", value = 10),
      numericInput("linkgenesize", "网络字体大小", value = 4),
      actionButton("submit", "提交")
    ),
    mainPanel(
      plotOutput("corrPlot")
    )
  )
)


server <- function(input, output) {
    observeEvent(input$submit, {
        req(input$datafile)

        # 读取数据文件
        #data <- read.delim(input$datafile$datapath, header = TRUE, check.names = FALSE,row.names = 1)
        tryCatch({
          data0 <- read.delim(input$datafile$datapath, header = TRUE, check.names = FALSE, row.names = 1)
        }, error = function(e) {
          print(e)
        })
        
        linkgene <- unlist(strsplit(input$linkgene, ",\\s*"))
        heatgene <- unlist(strsplit(input$heatgene, ",\\s*"))
        
        
        ge <- as.data.frame(data0[, linkgene, drop = FALSE])
        im <- as.data.frame(data0[, heatgene, drop = FALSE])
        write.table(colnames(im),"id.txt",sep="\t")
        rumen<-ge
        genus<-im

gene<-linkgene
exp<-t(data0)[union(gene,heatgene),]


dimnames=list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), ncol = ncol(exp), dimnames = list(rownames(exp), colnames(exp)))
data=avereps(data)
corFilter=0 #设置相关性系数阈值
pFilter=0.05#设置检验P值
namesize = 10#
linkgenesize<-4

corFilter <- input$corFilter
pFilter <- input$pFilter
namesize <- input$namesize
linkgenesize<-input$linkgenesize

#rt=as.data.frame(data[rowMeans(2^data-1)>1,])
rt<-data

#gene=c("HAT1","KAT2A") #设置需要寻求相关性基因的基因
#gene<- rownames(filter)

#xx <- apply(data[gene, , drop = FALSE], 2, as.numeric)

# 假设data[gene, , drop = FALSE]是转换前的数据
tempData <- data[gene, , drop = FALSE]

# 转换为数值型向量
xx <- apply(tempData, c(1,2), as.numeric)

# 重新赋予行名
rownames(xx) <- rownames(tempData)

rm(outTab)
rm(outTab2)
outTab2=data.frame()
outTab=data.frame()

for (n in 1:length(gene) ) {
  gene1=unlist(strsplit(gene,"\\|",))[n]
  x<-xx[gene1,]
  outTab=data.frame()
  for(j in rownames(rt)){
    y=as.numeric(rt[j,])
    #y=apply(log2(rt[j,]+1),2,as.numeric)
    #y=as.numeric(rt[j,])
    gene2=unlist(strsplit(j,"\\|",))[1]
    #gene2Type=unlist(strsplit(j,"\\|",))[2]
    #if(gene2Type=="protein_coding")
    corT=cor.test(x,y)
    gene1Name=unlist(strsplit(gene1,"\\|",))[1]
    gene2Name=unlist(strsplit(gene2,"\\|",))[1]
    
    z=lm(y~x)
    cor=corT$estimate
    cor=round(cor,3)
    pvalue=corT$p.value
    if(pvalue<0.05){
      pval=signif(pvalue,4)
      pval=format(pval,scientific=TRUE)
    }else{
      pval=round(pvalue,3)}
    
    if(abs(cor)>corFilter){
      outTab=rbind(outTab,cbind(gene1,gene2,cor,pvalue))
      outTab2=rbind(outTab2,outTab)
      outTab2<-outTab2[!duplicated(outTab2[c("gene1","gene2")]),]
    }}}

colnames(outTab2)<-c("spec","env","r","p")

write.table(outTab2,"cor—result.txt",sep="\t")
mantel<-outTab2
mantel$r<-as.numeric(mantel$r)

mantel <- mantel%>% 
  mutate(rd = cut(r, breaks = c(-Inf,0, Inf),
                  labels = c("< 0","> 0")),
rr=abs(r),
         pd = cut(rr, breaks = c(-Inf,0.3,0.5 ,Inf),
                  labels = c("0 - 0.3", "> 0.3", "> 0.5")))

#mantel <- mantel[mantel[,1] != "CR1", ]

  
output$corrPlot <- renderPlot({
  
qcorrplot(correlate(genus), type = "lower", diag = FALSE) +
  geom_square(color = NA) +  # 去除小格边框
  geom_mark(sep = '\n', size = 0, sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05, color = "white") +
  geom_couple(aes(colour = rd, size = pd), 
              data = mantel, 
              curvature = nice_curvature(),label.size = linkgenesize) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))) +  # 改变颜色方案
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#1B9E77","#D95F02")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(color = "black"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's r", 
                               override.aes = list(size = 1), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme(panel.border = element_blank(),  # 去除面板边框
        axis.line = element_blank(),    # 去除轴线
        text = element_text(size = 10),axis.text.x = element_text(size = namesize),
        axis.text.y = element_text(size = namesize))
    })
    })
}

shinyApp(ui = ui, server = server)

