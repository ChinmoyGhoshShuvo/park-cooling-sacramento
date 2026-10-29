
lib <- .libPaths()[1]
.libPaths(c(lib, .libPaths()))
pk <- c("terra","sf","MatchIt","cobalt","fixest","ggplot2","rstac","osmdata",
        "sensitivitymw","rbounds","data.table","jsonlite","httr2","exactextractr",
        "future.apply","patchwork","scales","marginaleffects")
need <- pk[!pk %in% rownames(installed.packages())]
cat("Installing:", need, "\n")
if (length(need)) install.packages(need, lib = lib, repos = "https://cloud.r-project.org",
                                   type = "binary", dependencies = c("Depends","Imports","LinkingTo"))
ok <- sapply(pk, requireNamespace, quietly = TRUE)
print(ok)
cat(if (all(ok)) "ALL_OK\n" else "MISSING\n")
