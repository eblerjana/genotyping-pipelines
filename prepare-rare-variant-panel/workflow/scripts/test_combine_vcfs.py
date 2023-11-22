import unittest, importlib
import tempfile
import io
import gzip
from contextlib import redirect_stdout
from combine_vcfs import VcfRecord, find_overlaps, create_combined_cluster


class Test_VcfRecord(unittest.TestCase):
	def test_vcfrecord(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"
		record = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		self.assertEqual(record.chrom(), "chr1")
		self.assertEqual(record.pos(), 1110696)
		self.assertEqual(record.id(), "rs6040355")
		self.assertEqual(record.ref(), "A")
		self.assertEqual(record.alt(), ["G", "T"])


	def test_same_variant(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"
		vcf_line2 = "chr1\t1110696\trs6040355\tA\tT,G\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"
		
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample4", "sample5", "sample6"])
		
		self.assertTrue(record1.same_variant_location(record2))
		self.assertTrue(record1.same_variant_location(record1))

		vcf_line3 = "chr1\t1110697\trs6040355\tA\tC\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t0|1:21:6:23,27\t1|1:2:0:18,2\t0/1:35:4"
		record3 = VcfRecord(vcf_line3, ["sample1", "sample2", "sample3"])
		self.assertFalse(record3.same_variant_location(record1))
		self.assertFalse(record3.same_variant_location(record2))

		vcf_line4 = "chr2\t1110697\trs6040355\tA\tC\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t0|1:21:6:23,27\t1|1:2:0:18,2\t0/1:35:4"
		record4 = VcfRecord(vcf_line4, ["sample1", "sample2", "sample3"])
		self.assertFalse(record4.same_variant_location(record3))
		self.assertFalse(record4.same_variant_location(record1))
		
	
	def test_combine_variants1(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"
		vcf_line2 = "chr1\t1110696\trs6040355\tA\tT,G\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t0|1:21:6:23,27\t1|1:2:0:18,2\t0/1:35:4"
		expected = "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT\t1|2\t2|1\t2/2\t0|2\t2|2\t0/2"
		
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample4", "sample5", "sample6"])
		expected_record = VcfRecord(expected, ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"])
		self.assertTrue(record1.same_variant_location(record2))

		record1.combine_variants(record2)
		self.assertTrue(expected_record == record1)
		self.assertTrue(expected_record.str_record() == record1.str_record())


	def test_combine_variants2(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"
		vcf_line2 = "chr2\t1110696\trs6040355\tA\tT\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t0|1:21:6:23,27\t1|1:2:0:18,2\t0/1:35:4"
		
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample4", "sample5", "sample6"])
		self.assertRaises(RuntimeError, record1.combine_variants, record2)


	def test_combine_variants3(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"
		vcf_line2 = "chr1\t1110696\trs6040355\tA\tT,G\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t0|1:21:6:23,27\t1|2:2:0:18,2\t0/2:35:4"
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample4", "sample5", "sample6"])

		record1.combine_variants(record2)
		expected = "chr1\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT\t1|2\t2|1\t2/2\t0|2\t2|1\t0/1"
		expected_record = VcfRecord(expected, ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"])

		self.assertTrue(expected_record == record1)


	def test_combine_variants4(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|2\t2|1\t2/2"
		vcf_line2 = "chr1\t1110696\trs6040355\tATTG\tT,GA\t67\tPASS\t.\tGT:GQ:DP:HQ\t0|1\t1|2\t0/2"
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample4", "sample5", "sample6"])

		record1.combine_variants(record2)
		expected = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT\t1|2\t2|1\t2/2\t0|2\t2|1\t0/1"
		expected_record = VcfRecord(expected, ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"])

		self.assertTrue(expected_record == record1)


	def test_combine_variants5(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|2\t2|1\t2/2"
		vcf_line2 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t0|1\t1|0\t0/0"
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample4", "sample5", "sample6"])

		record1.combine_variants(record2)
		expected = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT\t1|2\t2|1\t2/2\t0|1\t1|0\t0/0"
		expected_record = VcfRecord(expected,  ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"])
		self.assertTrue(expected_record == record1)


	def test_combine_variants6(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|2\t2|1\t2/2"
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		str_rec = record1.str_record()
		expected = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT\t1|2\t2|1\t2/2"
		self.assertTrue(str_rec == expected)

		str_rec = record1.str_record(sample_names = ["sample4", "sample5"])
		expected = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT\t1|2\t2|1\t2/2\t0/0\t0/0"
		self.assertTrue(str_rec == expected)
		
	def test_combine_variants7(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|2\t2|1\t2/2"
		vcf_line2 = "chr1\t1110696\trs6040355\tATTG\tCTG\t67\tPASS\t.\tGT:GQ:DP:HQ\t0|1\t1|0\t0/0"
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample2", "sample5", "sample6"])

		# should raise an error because samples are overlapping
		self.assertRaises(RuntimeError, record1.combine_variants, record2)


	def test_vcfrrecord_sort1(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|2\t2|1\t2/2"
		vcf_line2 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|0\t2|0\t1/2"
		vcf_line3 = "chr1\t1110696\trs6040355\tATTG\tA\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|2\t2|1\t2/2"
		
		variants = [VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"]), VcfRecord(vcf_line3, ["sample4", "sample5", "sample6"]), VcfRecord(vcf_line2, ["sample1", "sample2", "sample3"])]
		sorted_variants = sorted(variants)

		self.assertTrue(sorted_variants[0] == VcfRecord(vcf_line3, ["sample4", "sample5", "sample6"]))
		self.assertTrue(sorted_variants[1] == VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"]))
		self.assertTrue(sorted_variants[2] == VcfRecord(vcf_line2, ["sample1", "sample2", "sample3"]))


	def test_vcfrrecord_sort2(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|2\t2|1\t2/2"
		vcf_line2 = "chr1\t1110696\trs6040355\tATTG\tT,GA\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|0\t2|0\t1/2"

		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample1", "sample2", "sample3"])

		# both records are equal (in the sense that position + alleles are same)
		self.assertFalse(record1 < record2)
		self.assertFalse(record2 < record1)

		# record3 should be larger because the ALT allele is lexicographically larger
		vcf_line3 = "chr1\t1110696\trs6040355\tATTG\tT\t67\tPASS\t.\tGT:GQ:DP:HQ\t1|0\t2|0\t1/2"
		record3 = VcfRecord(vcf_line3, ["sample1", "sample2", "sample3"])
		self.assertTrue(record1 < record3)
		self.assertTrue(record2 < record3)


class Test_find_overlaps(unittest.TestCase):
	def test_find_overlaps1(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tGA,T\t67\tPASS\t.\tGT\t1|2\t2|1\t2/2"
		vcf_line2 = "chr1\t1110696\trs6040355\tATTG\tCTG\t67\tPASS\t.\tGT\t0|1\t1|0\t0/0"
		vcf_line3 = "chr1\t1110696\trs6040355\tATT\tCTG\t67\tPASS\t.\tGT\t0|1\t1|0\t0/0"
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample1", "sample2", "sample3"])
		record3 = VcfRecord(vcf_line3, ["sample3", "sample4", "sample5"])
		
		result = [r for r in find_overlaps([record1, record2, record3])]
		# no overlaps, since all variants are considered different
		self.assertTrue(len(result) == 3)
		
		result = [r for r in find_overlaps([record1, record1, record2])]
		self.assertTrue(len(result) == 2)
		self.assertTrue(result[0] == [record1, record1])
		self.assertTrue(result[1] == [record2])


class Test_create_combined_clusters(unittest.TestCase):
	def test_create_combine(self):
		vcf_line1 = "chr1\t1110696\trs6040355\tATTG\tCTG\t67\tPASS\t.\tGT\t1|1\t1|1\t1/1"
		vcf_line2 = "chr1\t1110696\trs6040355\tATTG\tCTG\t67\tPASS\t.\tGT\t0|1\t1|0\t0/0"
		vcf_line3 = "chr1\t1110696\trs6040355\tATT\tCTG\t67\tPASS\t.\tGT\t0|1\t1|0\t0/0"
		record1 = VcfRecord(vcf_line1, ["sample1", "sample2", "sample3"])
		record2 = VcfRecord(vcf_line2, ["sample4", "sample5", "sample6"])
		record3 = VcfRecord(vcf_line3, ["sample3", "sample4", "sample5"])

		clusters = [c for c in find_overlaps([record1, record2, record3])]
		result = []
		for c in clusters:
			result.append(create_combined_cluster(c))


		self.assertTrue(len(result) == 2)
		expected_line1 = "chr1\t1110696\trs6040355\tATTG\tCTG\t67\tPASS\t.\tGT\t1|1\t1|1\t1/1\t0|1\t1|0\t0/0"
		expected_line2 = "chr1\t1110696\trs6040355\tATT\tCTG\t67\tPASS\t.\tGT\t0/0\t0/0\t0|1\t1|0\t0/0\t0/0"
		computed_line1 = result[0].str_record(["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"])
		computed_line2 = result[1].str_record(["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"])

		self.assertTrue(expected_line1 == computed_line1)
		self.assertTrue(expected_line2 == computed_line2)


if __name__ == '__main__':
	unittest.main()
