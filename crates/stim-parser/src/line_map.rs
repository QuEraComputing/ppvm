/// Maps byte offsets in source to 1-indexed line/column positions.
pub struct LineMap {
    /// `starts[i]` = byte offset of the start of line (i+1).
    starts: Vec<usize>,
}

impl std::fmt::Debug for LineMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LineMap")
            .field("lines", &self.starts.len())
            .finish()
    }
}

impl LineMap {
    /// Build a `LineMap` for `src`.
    pub fn new(src: &str) -> Self {
        let mut starts = vec![0];
        for (i, b) in src.bytes().enumerate() {
            if b == b'\n' {
                starts.push(i + 1);
            }
        }
        Self { starts }
    }

    /// 1-indexed line number for a byte offset.
    pub fn line_of(&self, byte_offset: usize) -> usize {
        match self.starts.binary_search(&byte_offset) {
            Ok(i) => i + 1,
            Err(i) => i, // i is the insertion index; start of line `i` is at starts[i-1].
        }
    }

    /// 1-indexed `(line, col)` for a byte offset.
    pub fn line_col(&self, byte_offset: usize) -> (usize, usize) {
        let line = self.line_of(byte_offset);
        let line_start = self.starts[line - 1];
        let col = byte_offset - line_start + 1;
        (line, col)
    }

    /// Byte offset of the start of line `(line_idx + 1)`. `None` for
    /// out-of-range indices.
    pub fn starts_at(&self, line_idx: usize) -> Option<usize> {
        self.starts.get(line_idx).copied()
    }
}

#[cfg(test)]
mod line_map_tests {
    use super::LineMap;

    #[test]
    fn line_col_at_line_start() {
        let m = LineMap::new("abc\ndef\nghi");
        assert_eq!(m.line_col(0), (1, 1));
        assert_eq!(m.line_col(4), (2, 1));
        assert_eq!(m.line_col(8), (3, 1));
    }

    #[test]
    fn line_col_mid_line() {
        let m = LineMap::new("abc\ndef\nghi");
        assert_eq!(m.line_col(2), (1, 3));
        assert_eq!(m.line_col(6), (2, 3));
    }

    #[test]
    fn line_col_at_eof() {
        let m = LineMap::new("abc\ndef");
        assert_eq!(m.line_col(7), (2, 4));
    }
}
