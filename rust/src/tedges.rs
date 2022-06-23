#[derive(Clone,Copy,Debug, Eq)]
pub struct TEdge (pub u8, pub u8, pub u8);
impl PartialEq for TEdge {
    fn eq(&self, other: &Self) -> bool {
        (self.2 == other.2) &&
            ((self.0 == other.0 && self.1 == other.1) ||
                self.0 == other.1 && self.1 == other.0)
    }
}

