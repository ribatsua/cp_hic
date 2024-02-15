htmlheader <- function() {
  cat("<!DOCTYPE html>\n")
  cat("<html lang=\"en\">\n")
  cat("<head>\n")
  cat("    <meta charset=\"utf-8\">\n")
  cat("    <meta content=\"IE=edge\" http-equiv=\"X-UA-Compatible\">\n")
  cat("    <meta content=\"width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no\" name=\"viewport\">\n")
  cat("    <meta content=\"\" name=\"description\">\n")
  cat("    <meta content=\"\" name=\"", paste0(Sys.info()["user"][[1]]), "\">\n")
  cat("    <link href=https://igv.org/web/release/2.9.1/examples/img/favicon.ico rel=\"igv icon\">\n")
  cat("    <title>IGV - Dev</title>\n")
  cat("    <script src=\"https://cdn.jsdelivr.net/npm/igv@2.9.1/dist/igv.min.js\"></script>\n")
  cat("</head>\n")
  
  cat("<body>\n")
  cat("<h1><u>higashi compare dcHic from Hi-C data</u></h1>\n")
  cat("<p style=\"color: #C86400\"><b>A compartment</b></p>\n")
  cat("<p style=\"color: #0064C8\"><b>B compartment</b></p>\n")
  cat("<p style=\"color: #99FFCC\"><b>Mahalanobis distance: Mahalanobis distance score to represent outlierness of the bin</b></p>\n")
  cat("<p style=\"color: #009900\"><b>log10Padj: -log10 of the P.adjusted value of corresponding Mdist score</b></p>\n")
  cat("<div id=\"igvDiv\" style=\"padding-top: 10px;padding-bottom: 10px; border:1px solid lightgray\"></div>\n")
}
scriptbody <- function(genome,start=F,end=F) {
  
  if (start) {
    cat ("<script type=\"text/javascript\">
    document.addEventListener(\"DOMContentLoaded\", function () {
      var igvDiv = document.getElementById(\"igvDiv\");
      var options = {
        locus: '19:49301000-49305700',
        genome:",paste0("\"",genome,"\""),",
	tracks: [\n")
  }
  if (end) {
    cat ("
       ]
     }; igv.createBrowser(igvDiv, options)
                .then(function (browser) {
                    console.log(\"Created IGV browser\")});
   })
  </script>
 </body>
</html>\n")
  }
}

htmlbody <- function(bedgraph_paths, genome) {
  
  scriptbody(genome, start = TRUE)
  
  for (file_path in bedgraph_paths) {

    
    
    # Process the rest of the file and generate tracks...
    
    # Output BedGraph track configuration in HTML
    #cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
#track name=\"", colnames(compbdg)[j], " PC\" description=\"BedGraph format\" visibility=full color=", color, " altColor=", altcolor, " priority=20 plotType=\"points\"\n",
        #file = paste0(file_path, ".PC.bedGraph"))
    
    #write.table(compbdg[, c(1:3, j)], file = paste0(file_path, ".PC.bedGraph"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    #R.utils::gzip(paste0(file_path, ".bedGraph"), paste0(file_path, ".PC.bedGraph.gz"), overwrite = TRUE)
    # R.utils::gzip(file_path, paste0(file_path, ".PC.bedGraph.gz"), overwrite = TRUE)
    # bed_name<- regmatches(file_path, regexpr("(?<=\\/)[^\\/]+(?=\\.)", file_path, perl = TRUE))
    bed_name<-regmatches(file_path, regexpr("(?<=\\/)[^\\/]+(?=\\.[^.]*$)", file_path, perl = TRUE))
    # bed_name<- sub(".*/([^/.]+)\\..*", "\\1", file_path)
    if (!endsWith(file_path, ".gz")) {
    # 如果不是，使用 R.utils::gzip 进行压缩
    R.utils::gzip(file_path, paste0(file_path, ".gz"), overwrite = TRUE)
    file_path <- paste0(file_path, ".gz")
  }
    # cmd <- paste0("create_datauri ",file_path,".PC.bedGraph.gz")
    cmd <- paste0("create_datauri ",file_path)
    datauri<- file_path
    datauri <- system(cmd, wait = TRUE, intern = TRUE)
    #datauri<-paste0(file_path, ".PC.bedGraph.gz")
    cat("
          {
              name: ", paste0("\'", bed_name, ".PC',"), "
              url: ",paste0("\"",datauri,"\""),",
              indexed: false,
              format: \"bedGraph\"
          },\n")
  }
  scriptbody(genome, end = TRUE)
}


html_file_path <- paste0("/home/dchic/dcHiC/snm3c250k/DifferentialResult/higashi_ori_prebulk_compare_dchic.html")

folder_path<-'/home/dchic/dcHiC/snm3c250k/DifferentialResult/snm3c/viz/vizIGV_intra/data'

# 获取文件夹下所有文件的路径
bedGraph_file_path <- list.files(folder_path, full.names = TRUE)


genome<-'hg19'


sink(html_file_path)
# 生成 HTML 头部
htmlheader()

# 生成 HTML 主体
htmlbody(bedGraph_file_path, genome = 'hg19')

# 恢复默认的输出
sink()

cat("HTML file has been saved to:", html_file_path, "\n")