test_check_parser <-function()
{
    b <- PGA:::check_parser()
    checkEquals(b,TRUE)
}
