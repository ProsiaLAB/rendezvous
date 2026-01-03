pub fn murmur3_32(key: &[u8], seed: u32) -> u32 {
    let c1: u32 = 0xcc9e2d51;
    let c2: u32 = 0x1b873593;
    let mut h1: u32 = seed;
    let len = key.len() as u32;
    let nblocks = len / 4;

    // body
    for i in 0..nblocks {
        let mut k1: u32 = u32::from_le_bytes([
            key[(i * 4) as usize],
            key[(i * 4 + 1) as usize],
            key[(i * 4 + 2) as usize],
            key[(i * 4 + 3) as usize],
        ]);

        k1 = k1.wrapping_mul(c1);
        k1 = k1.rotate_left(15);
        k1 = k1.wrapping_mul(c2);

        h1 ^= k1;
        h1 = h1.rotate_left(13);
        h1 = h1.wrapping_mul(5).wrapping_add(0xe6546b64);
    }

    // tail
    let mut k1: u32 = 0;
    let tail_index = (nblocks * 4) as usize;
    let tail_size = (len % 4) as usize;

    if tail_size >= 3 {
        k1 ^= (key[tail_index + 2] as u32) << 16;
    }
    if tail_size >= 2 {
        k1 ^= (key[tail_index + 1] as u32) << 8;
    }
    if tail_size >= 1 {
        k1 ^= key[tail_index] as u32;
        k1 = k1.wrapping_mul(c1);
        k1 = k1.rotate_left(15);
        k1 = k1.wrapping_mul(c2);
        h1 ^= k1;
    }

    // finalization
    h1 ^= len;
    h1 ^= h1 >> 16;
    h1 = h1.wrapping_mul(0x85ebca6b);
    h1 ^= h1 >> 13;
    h1 = h1.wrapping_mul(0xc2b2ae35);
    h1 ^= h1 >> 16;

    h1
}
