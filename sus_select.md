# Introduction #

We need a solution s.t. SUS selection can be done in constant time...

i.e. `chrom_t sus_select(float target_f)` is constant

Given `chrom_t` ordered in `vector<chrom_t>` with `c` members,
we want to return
`vector[i] s.t. cum_f[i] <= target_f < cum_f[i+1]`
where `cum_t[i] = sum{n=0,n=i)}(vector[n])`

`cum_f[i]` can be precomputed by:
`cum_f[0] = vector[0].f();`
`cum_f[n] = cum_f[n-1] + vector[n].f();`

We want O(1) algorithm that will return `i` in `cum_f` as defined above.

O(log c) algorithm is binary search over `cum_f`

`:------:---:--:--------------:----:-:--:`

If numbers were integers, could create a memory block of
`cum_f[c-1]` integers, called `lookup`. `lookup[0]->lookup[cum_f[0]]`
would contain the value 0. `lookup[cum_f[0]]->lookup[cum_f[1]]`
would contain the value 1. Generally:
`lookup[cum_f[n-1]]->lookup[cum_f[n]]` would contain the value `n`.
Hence `i = lookup[target_f]`

There are RLE encodings that allow searching. Ask Pete about this...

Problem: We are not using integers :-(

Can convert to integers by finding Lowest Common Denominator (LCD) of
fitnesses, and assigning this value to 1. All fitnesses can then be scaled
to a multiple of this.

Problem: LCD compilation is slow.
Problem: LCD of 1/999 and 1/1000 is 1/999000 which will lead
to large amounts of memory storage.
Problem: Floating-point is not exact, especially when finding
LCD. Therefore cannot guarantee all results are exact multiples
of LCD. Therefore cannot guarantee all results are integer
without rounding. We have to accept rounding in floating point
anyway.

```

:------:---:--:--------------:----:-:--:

:-----
-:---:
--:---
------
-----:
----:-
:--:

0
1
2
3
4
5
6
```
```
row = floor(target_f)
if(row.start)==(row.end)
  i = row.start
else
  i = row.search(target_f);
end
return cum_f[i];
```
To construct rows:
(And we can pre-scale the data to make a good fit)
```
n = 0;
for(i = 0; i < floor(cum_f[c-1]); i++)
  if(n==c || cum_f[n+1] > i)
    row.start = row.end = n;
  else
    row.start = n;
    for(; (n <= c) && (cum_f[n] < i); ++n)
      row.add(cum_f[n],n);
    row.end = n;
  end
end
```
where
```
class row {
  public:
  int start, end;
  void add(float);
  int search(float);
  
  private:
  map<float, int> data;
};
```