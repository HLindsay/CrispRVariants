context("Counting mutations")

alns <- GenomicAlignments::GAlignments(seqnames=rep("1", 5), 
                                       pos=as.integer(rep(1,5)), 
                    cigar=c("10M", "5M1D5M", "4M1D6M4D2M",
                            "1M2D3M4I", "2M3I2M3I2M"), 
                    strand = S4Vectors::Rle(factor(rep("+", 5),
                                         levels = c("+","-","*"))))

test_that("countDeletions returns the expected counts", {
  expect_equal(countDeletions(alns), 1)
  expect_equal(countDeletions(alns, multi.del = TRUE), 2)
  expect_equal(countDeletions(alns, del.and.ins = TRUE), 2)
  expect_equal(countDeletions(alns, multi.del = TRUE, del.and.ins = TRUE), 3) 
  expect_equal(countDeletions(alns, del.ops=c("N")), 0) 
})

test_that("countInsertions returns the expected counts", {
  expect_equal(countInsertions(alns), 0)
  expect_equal(countInsertions(alns, multi.ins = TRUE), 1)
  expect_equal(countInsertions(alns, ins.and.del = TRUE), 1)
  expect_equal(countInsertions(alns, multi.ins = TRUE, ins.and.del = TRUE), 2) 
})

test_that("countIndels returns the expected counts", {
  expect_equal(countIndels(alns), 4)
})

test_that("indelPercent returns the expected value", {
  expect_equal(indelPercent(alns), (4/5)*100)
})


context("Accessory plotting functions")

test_that("Extraploation of x-tick labels works correctly", {

    tck <- .getAxisCoords(c(1:20), label.at = 10, lab.boundaries = c(-1,1))
    # Only boundaries within the labels are included 
    expect_equal(tck$tick_labs, c(1,10,20))
    
    # Error if labels aren't the same length as locations
    expect_error(.getAxisCoords(c(1:20), "A"))  
    
    tck <- .getAxisCoords(c(1:23), labels = c(-17:-1, 1:6),
                          label.at = 2, lab.boundaries = c(-1,1))
    
    tmp <- setdiff(c(seq(-16,6, by = 2), c(-1,1)),0)
    exp_result <- tmp[order(tmp)] 
    expect_equal(tck$tick_labs, exp_result)

    # If min.tick.sep is increased, ticks close to boundary
    # should be removed
    tck <- .getAxisCoords(c(1:23), labels = c(-17:-1, 1:6),
                          label.at = 2, lab.boundaries = c(-1,1),
                          min.tick.sep = 2)
    expect_equal(tck$tick_labs, exp_result[!exp_result %in% c(-2,2)])
})
