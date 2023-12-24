package com.hartwig.hmftools.esvee.models;

import org.jetbrains.annotations.Nullable;

public interface TrimmableAssembly<SELF>
{
    @Nullable
    SELF trim(final int removeLeft, final int removeRight);
}