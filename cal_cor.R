





#' 对数据框进行相关性计算
#'
#' @param df_x 数据框x
#' @param df_y 数据框y，为NULL则数据框x自己跟自己计算相关性。存在则数据框x与y计算相关性
#' @param cor_type 可选择 "pearson" "spearman" "lm"
#' @param alternative
#'
#'
#' @importFrom rstatix cor_gather
#' @importFrom broom glance
#' @import purrr
#' @import dplyr
#'
#'
#' @return 返回相关系数和pvalue
#' @export cal_cor
#'
#' @examples # cal_cor(df_x = mtcars,df_y = mtcars)
#'
#'
#'
cal_cor <- function(df_x,df_y = NULL,cor_type = "pearson",alternative = "two.sided"){

  # library(rstatix) cor_gather
  # library(broom)   glance
  library(purrr)
  library(dplyr)

  # 有数据框x和数据框y的情况
  if (is.null(df_y)) {
    df_y <- df_x
  }


  if (cor_type %in% c("pearson","spearman")){

    # 计算相关系数
    r <- cor(x = df_x,y = df_y,method = cor_type)

    # 计算相关系数对应的p值
    p <- map_df(df_x ,function(x){
      map(df_y,function(y){
        cor.test(as.numeric(x),method = cor_type,
                 as.numeric(y))$p.value
      }
      )
    }
    ) %>% as.data.frame() %>% `rownames<-`(colnames(df_x))

    # 整合pvalue和cor 并计算adj.p
    rp_long_df <- inner_join(rstatix::cor_gather(r),rstatix::cor_gather(p) %>% dplyr::rename(p = cor),by = c("var1",  "var2")) %>%
      mutate(p.adj = p.adjust(p,method = "BH"))
  }


  if (cor_type == "lm") {

    # 计算拟合R方
    r <- map_df(df_x ,function(x){
      map(df_y,function(y){
        data = data.frame(x = as.numeric(x),y = as.numeric(y))
        model = lm(x~y, data = data)
        r.squared = summary(model)$r.squared
        r.squared
      }
      )
    }
    ) %>% as.data.frame() %>% `rownames<-`(colnames(df_x))

    #计算拟合R方的pvalue
    p <- map_df(df_x ,function(x){
      map(df_y,function(y){
        data = data.frame(x = as.numeric(x),y = as.numeric(y))
        model = lm(x~y, data = data)
        p.value = broom::glance(model)$p.value
        p.value
      }
      )
    }
    ) %>% as.data.frame() %>% `rownames<-`(colnames(df_x))


    rp_long_df <- inner_join(rstatix::cor_gather(r),rstatix::cor_gather(p) %>% dplyr::rename(p = cor),by = c("var1",  "var2")) %>%
      mutate(p.adj = p.adjust(p,method = "BH")) %>%
      dplyr::rename(R.square = cor)

  }


  out = list(rp = rp_long_df,r = r,p = p)
  return(out)

}
