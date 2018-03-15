BEGIN_PROVIDER [ character*(32), DFT_TYPE]
 implicit none
 BEGIN_DOC
! defines the type of DFT applied: LDA, GGA etc ... 
 END_DOC
 logical :: is_short_range_lda,is_lda
 if(correlation_functional.eq."None")then
  is_short_range_lda = (exchange_functional.EQ."short_range_LDA")
  is_lda             = (exchange_functional.EQ."LDA")
 else if(exchange_functional.eq."None")then
  is_short_range_lda = (correlation_functional.EQ."short_range_LDA")
  is_lda             = (correlation_functional.EQ."LDA")
 else 
  is_short_range_lda = (exchange_functional.EQ."short_range_LDA").and.(correlation_functional.EQ."short_range_LDA")
  is_lda = (exchange_functional.EQ."LDA").and.(correlation_functional.EQ."LDA")
 endif
 if(is_short_range_lda.or.is_lda)then
  DFT_TYPE = "LDA"
 else
  DFT_TYPE = "GGA"
 endif
END_PROVIDER 
